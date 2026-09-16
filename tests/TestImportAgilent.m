classdef TestImportAgilent < matlab.unittest.TestCase
    %TESTIMPORTAGILENT Unit tests for the standalone Agilent importers
    %   importAgilentUV, importAgilentCH, importAgilentMS.
    %
    % The tests run over every Agilent ".D" data folder found (recursively)
    % under a dataset directory. They fall into two groups:
    %
    %   * SELF-CONTAINED checks that need no external reference and therefore
    %     run on ANY .D folder: each file imports without error, all structs
    %     from one importer share a fixed field layout (so they concatenate
    %     into an array), field types are stable, the .CH header min/max match
    %     the decoded signal, and scaling is consistent.
    %
    %   * REFERENCE checks that compare the decoded data bit-for-bit against
    %     the CRAN chromConverter R package. These are skipped (marked
    %     "incomplete") unless the reference CSVs produced by
    %     tests/generate_groundtruth.R are present.
    %
    % Configuration (see README):
    %   * Dataset folder: environment variable AGILENT_TEST_DATA, else the
    %     "datasets" folder next to the repo root.
    %   * Reference folder: environment variable AGILENT_TEST_TRUTH, else
    %     <tempdir>/ca_agilent_truth (the default output of the R script).
    %
    % Run with:  runtests("tests")   or   runtests("tests/TestImportAgilent.m")

    properties
        DatasetDir (1,1) string = ""
        TruthDir   (1,1) string = ""
        DFolders   (1,:) string = string.empty
    end

    methods (TestClassSetup)
        function setUpPathsAndData(testCase)
            import matlab.unittest.fixtures.PathFixture

            here = fileparts(mfilename("fullpath"));
            repoRoot = fileparts(here);
            % Put the importers (repo root) on the path for the test's lifetime.
            testCase.applyFixture(PathFixture(repoRoot));

            dataDir = string(getenv("AGILENT_TEST_DATA"));
            if dataDir == ""
                dataDir = fullfile(repoRoot, "datasets");
            end
            testCase.assumeTrue(isfolder(dataDir), ...
                "Dataset folder not found: " + dataDir + newline + ...
                "Set the AGILENT_TEST_DATA environment variable to a folder " + ...
                "containing Agilent .D directories (see README).");
            testCase.DatasetDir = dataDir;

            d = dir(fullfile(dataDir, "**", "*.D"));
            d = d([d.isdir]);
            testCase.assumeNotEmpty(d, "No .D folders found under " + dataDir);
            testCase.DFolders = string(fullfile({d.folder}, {d.name}));

            td = string(getenv("AGILENT_TEST_TRUTH"));
            if td == ""
                td = fullfile(tempdir, "ca_agilent_truth");
            end
            testCase.TruthDir = td;
        end
    end

    methods (Test)   % ---------- self-contained checks ----------

        function uvImportsSchemaTypes(testCase)
            files = testCase.filesWithExt("uv");
            testCase.assumeNotEmpty(files, "No .UV files found.");
            structs = cell(1, numel(files));
            for i = 1:numel(files)
                raw = importAgilentUV(files(i).path, ApplyScaling=false);
                sc  = importAgilentUV(files(i).path);
                structs{i} = sc;
                testCase.verifyClass(sc.method.datetime, "datetime");
                testCase.verifyClass(sc.sample.vialPosition, "string");
                testCase.verifyGreaterThan(numel(sc.signal.wavelength), 0);
                testCase.verifyEqual(sc.signal.absorbance, ...
                    raw.signal.absorbance * raw.signal.scalingFactor, ...
                    "Scaled absorbance is not raw * scalingFactor.");
            end
            testCase.verifyFixedSchema(structs, "UV");
        end

        function chImportsSchemaTypesSelfCheck(testCase)
            files = testCase.filesWithExt("ch");
            testCase.assumeNotEmpty(files, "No .CH files found.");
            structs = cell(1, numel(files));
            for i = 1:numel(files)
                raw = importAgilentCH(files(i).path, ApplyScaling=false);
                sc  = importAgilentCH(files(i).path);
                structs{i} = sc;
                testCase.verifyClass(sc.method.datetime, "datetime");
                % Header-recorded min/max must match the decoded signal (v130).
                if ~isnan(raw.signal.minIntensity)
                    testCase.verifyEqual(min(raw.signal.intensity), raw.signal.minIntensity, ...
                        "Decoded minimum does not match the header for " + files(i).path);
                    testCase.verifyEqual(max(raw.signal.intensity), raw.signal.maxIntensity, ...
                        "Decoded maximum does not match the header for " + files(i).path);
                end
                % For v130 the intercept is 0, so scaled == raw * factor exactly.
                if raw.file.version == "130"
                    testCase.verifyEqual(sc.signal.intensity, ...
                        raw.signal.intensity * raw.signal.scalingFactor, ...
                        "Scaled intensity is not raw * scalingFactor for v130.");
                end
            end
            testCase.verifyFixedSchema(structs, "CH");
        end

        function msImportsSchemaTypes(testCase)
            files = testCase.filesWithExt("MS");
            testCase.assumeNotEmpty(files, "No .MS files found.");
            structs = cell(1, numel(files));
            for i = 1:numel(files)
                ms = importAgilentMS(files(i).path);
                structs{i} = ms;
                testCase.verifyClass(ms.method.datetime, "datetime");
                testCase.verifyEqual(numel(ms.signal.tic), numel(ms.signal.retentionTime));
                testCase.verifyEqual(size(ms.signal.xic), ...
                    [numel(ms.signal.retentionTime), numel(ms.signal.mz)]);
                testCase.verifyTrue(issorted(ms.signal.mz) && ...
                    isequal(ms.signal.mz, unique(ms.signal.mz)), ...
                    "m/z axis is not sorted-unique.");
            end
            testCase.verifyFixedSchema(structs, "MS");
        end
    end

    methods (Test)   % ---------- reference checks (chromConverter) ----------

        function uvMatchesChromConverter(testCase)
            testCase.assumeReferenceAvailable();
            files = testCase.filesWithExt("uv");
            matched = false;
            for i = 1:numel(files)
                tp = fullfile(testCase.TruthDir, "uv_" + files(i).key + ".csv");
                if ~isfile(tp), continue; end
                matched = true;
                s = importAgilentUV(files(i).path, ApplyScaling=false);
                T = readmatrix(tp, NumHeaderLines=1);                 % col 1 = rt, rest = wavelengths
                testCase.verifyEqual(size(s.signal.absorbance), size(T(:, 2:end)), ...
                    "UV matrix size differs from chromConverter for " + files(i).path);
                testCase.verifyEqual(max(abs(s.signal.absorbance - T(:, 2:end)), [], "all"), 0, ...
                    "UV data differs from chromConverter for " + files(i).path);
                testCase.verifyLessThan(max(abs(s.signal.retentionTime - T(:, 1))), 1e-9);
            end
            testCase.assumeTrue(matched, "No UV reference CSVs matched the discovered files.");
        end

        function chMatchesChromConverter(testCase)
            testCase.assumeReferenceAvailable();
            files = testCase.filesWithExt("ch");
            matched = false;
            for i = 1:numel(files)
                tp = fullfile(testCase.TruthDir, "ch_" + files(i).key + ".csv");
                if ~isfile(tp), continue; end
                matched = true;
                s = importAgilentCH(files(i).path, ApplyScaling=false);
                T = readmatrix(tp, NumHeaderLines=1);                 % rt, raw
                testCase.verifyEqual(numel(s.signal.intensity), size(T, 1), ...
                    "CH length differs from chromConverter for " + files(i).path);
                testCase.verifyEqual(max(abs(s.signal.intensity - T(:, 2))), 0, ...
                    "CH data differs from chromConverter for " + files(i).path);
            end
            testCase.assumeTrue(matched, "No CH reference CSVs matched the discovered files.");
        end

        function msMatchesChromConverter(testCase)
            testCase.assumeReferenceAvailable();
            files = testCase.filesWithExt("MS");
            summaryPath = fullfile(testCase.TruthDir, "ms_summary.csv");
            haveSummary = isfile(summaryPath);
            if haveSummary
                summary = readtable(summaryPath, TextType="string");
            end
            matched = false;
            for i = 1:numel(files)
                tp = fullfile(testCase.TruthDir, "ms_" + files(i).key + ".csv");
                if ~isfile(tp), continue; end
                matched = true;
                ms = importAgilentMS(files(i).path);
                T = readmatrix(tp, NumHeaderLines=1);                 % rt, tic
                testCase.verifyEqual(numel(ms.signal.tic), size(T, 1), ...
                    "MS scan count differs from chromConverter for " + files(i).path);
                testCase.verifyLessThan(max(abs(ms.signal.retentionTime - T(:, 1))), 1e-9);
                testCase.verifyEqual(max(abs(ms.signal.tic - T(:, 2))), 0, ...
                    "MS TIC differs from chromConverter for " + files(i).path);
                if haveSummary
                    row = summary(summary.key == files(i).key, :);
                    if height(row) == 1
                        testCase.verifyEqual(sum(ms.signal.xic, "all"), row.intensity_sum, ...
                            "MS abundance sum differs from chromConverter for " + files(i).path);
                        testCase.verifyEqual(max(ms.signal.xic, [], "all"), row.intensity_max, ...
                            "MS abundance max differs from chromConverter for " + files(i).path);
                    end
                end
            end
            testCase.assumeTrue(matched, "No MS reference CSVs matched the discovered files.");
        end
    end

    methods (Test)   % ---------- version-coverage smoke tests ----------

        function chAllVersionsRunAndShareSchema(testCase)
            % Relabel a real .CH file as each supported version. The data is
            % nonsense under the wrong decoder, so this only asserts that every
            % version code path runs and returns the same fixed schema.
            files = testCase.filesWithExt("ch");
            testCase.assumeNotEmpty(files, "No .CH files found.");
            ref = "";
            for ver = ["8" "30" "81" "130" "179" "181"]
                relabelled = TestImportAgilent.relabelVersion(files(1).path, ver);
                w = warning("off", "importAgilentCH:minMaxMismatch");
                cleanup = onCleanup(@() warning(w));
                s = importAgilentCH(relabelled);
                clear cleanup
                sig = TestImportAgilent.schemaSig(s);
                if ref == "", ref = sig; end
                testCase.verifyEqual(sig, ref, "CH version " + ver + " returns a different schema.");
            end
        end

        function msGcVariantDecodes(testCase)
            % Relabel a real LC-MS file as GC-MS (with the scan count also
            % written little-endian at 0x142) and check it is detected and
            % decodes the same data body.
            files = testCase.filesWithExt("MS");
            testCase.assumeNotEmpty(files, "No .MS files found.");
            [gcPath, lcPath] = TestImportAgilent.relabelMsGc(files(1).path);
            gc = importAgilentMS(gcPath);
            lc = importAgilentMS(lcPath);
            testCase.verifyEqual(gc.file.type, "GC / MS Data File");
            testCase.verifyEqual(gc.signal.xic, lc.signal.xic, ...
                "GC-relabelled file decodes a different body than the LC original.");
        end
    end

    methods (Access = private)
        function files = filesWithExt(testCase, ext)
            % Struct array of every file with the given extension across all .D
            % folders, with fields: path, dfolder, key (matching the R script).
            files = struct("path", {}, "dfolder", {}, "key", {});
            for df = testCase.DFolders
                listing = dir(fullfile(df, "*." + ext));
                for j = 1:numel(listing)
                    p = fullfile(listing(j).folder, listing(j).name);
                    files(end+1) = struct("path", string(p), "dfolder", df, ...
                        "key", TestImportAgilent.keyFor(df, p)); %#ok<AGROW>
                end
            end
        end

        function assumeReferenceAvailable(testCase)
            testCase.assumeTrue(isfolder(testCase.TruthDir), ...
                "Reference folder not found: " + testCase.TruthDir + newline + ...
                "Generate it with tests/generate_groundtruth.R (see README), " + ...
                "or set AGILENT_TEST_TRUTH.");
        end

        function verifyFixedSchema(testCase, structs, label)
            sigs = cellfun(@TestImportAgilent.schemaSig, structs);
            testCase.verifyTrue(isscalar(unique(sigs)), ...
                label + " schema differs across files (structs would not concatenate).");
            testCase.verifyTrue(TestImportAgilent.canConcat(structs), ...
                label + " structs do not concatenate into an array.");
        end
    end

    methods (Static, Access = private)
        function k = keyFor(dfolder, filePath)
            % Must match generate_groundtruth.R: sanitized ".D name" __ "file name".
            san = @(s) regexprep(string(s), "[^A-Za-z0-9]", "_");
            [~, dn, de] = fileparts(string(dfolder));
            [~, fn, fe] = fileparts(string(filePath));
            k = san(dn + de) + "__" + san(fn + fe);
        end

        function s = schemaSig(d)
            groups = string(fieldnames(d))';
            sub = arrayfun(@(g) strjoin(string(fieldnames(d.(g)))', ","), groups);
            s = strjoin([groups, sub], "|");
        end

        function tf = canConcat(structCell)
            try
                arr = structCell{1};
                for i = 2:numel(structCell), arr(i) = structCell{i}; end
                tf = numel(arr) == numel(structCell);
            catch
                tf = false;
            end
        end

        function out = relabelVersion(srcPath, ver)
            raw = TestImportAgilent.readBytes(srcPath);
            v = char(ver);
            raw(1) = numel(v);
            raw(2:1+numel(v)) = uint8(v);
            out = fullfile(tempdir, "ca_ch_v" + ver + ".ch");
            TestImportAgilent.writeBytes(out, raw);
        end

        function [gcPath, lcPath] = relabelMsGc(srcPath)
            raw = TestImportAgilent.readBytes(srcPath);
            % LC scan count: uint16 big-endian at 0x118 (bytes 281..282, 1-based).
            nScans = double(raw(281))*256 + double(raw(282));
            % Write it little-endian at 0x142 (bytes 323 low, 324 high).
            raw(323) = uint8(mod(nScans, 256));
            raw(324) = uint8(floor(nScans/256));
            label = uint8('GC / MS Data File');
            raw(5) = numel(label);
            raw(6:5+numel(label)) = label;
            gcPath = fullfile(tempdir, "ca_ms_GC.MS");
            TestImportAgilent.writeBytes(gcPath, raw);
            lcPath = fullfile(tempdir, "ca_ms_LCref.MS");
            copyfile(srcPath, lcPath);
        end

        function raw = readBytes(p)
            fid = fopen(p, "r"); raw = fread(fid, Inf, "uint8=>uint8"); fclose(fid);
        end

        function writeBytes(p, raw)
            fid = fopen(p, "w"); fwrite(fid, raw, "uint8"); fclose(fid);
        end
    end
end
