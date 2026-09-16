function fileContent = importAgilentCH(filePath, options)
%IMPORTAGILENTCH Import an Agilent ChemStation/OpenLab .CH signal file.
%
% Syntax:
%   fileContent = importAgilentCH(filePath)
%   fileContent = importAgilentCH(filePath, ApplyScaling=false)
%
% Description:
%   Reads a single-channel chromatogram (e.g. DAD1A.ch) from an Agilent .D
%   data folder: one intensity value per retention time.
%
%   A .CH file is NOT a spectrum. It is one signal channel the detector
%   produced. For a diode-array detector the channel is an absorbance averaged
%   over a wavelength band and corrected against a reference band; the exact
%   definition is stored in the file and returned as
%   'instrument.signalDescription', e.g. "DAD1A, Sig=254,16  Ref=450,100"
%   (254 nm over a 16 nm band, referenced against 450 nm / 100 nm band). Other
%   detectors (MWD, ELSD, FID) store their own signal in the same container.
%
%   For a diode-array run this file is complementary to, not redundant with,
%   the .UV spectral cube: the .CH channels typically cover the full run while
%   the .UV cube may be stored only over a shorter window, so the .CH traces
%   are not reconstructible from the .UV alone.
%
%   Supported file versions
%   -----------------------
%   Detected from the version tag at the start of the file:
%     "130" : ChemStation/OpenLab DAD/MWD/ELSD, delta-encoded      [VALIDATED]
%     "30"  : legacy ChemStation DAD/MWD/ELSD, delta-encoded       [UNTESTED]
%     "8"   : legacy FID/ADC, delta-encoded                        [UNTESTED]
%     "81"  : FID/ADC, double-delta-encoded                        [UNTESTED]
%     "181" : FID/ADC, double-delta-encoded                        [UNTESTED]
%     "179" : FID/ADC, uncompressed little-endian doubles          [UNTESTED]
%   Only version 130 files were available to test. The other five are
%   implemented from the chromConverter R package and rainbow project format
%   references and are marked UNTESTED; their decoders (double-delta and
%   double-array) and header offsets were ported but not verified against real
%   files. An unrecognised version tag is a hard error.
%
%   Version 179 has two sub-forms selected from the file-type code at offset
%   347 and, for GC files, the software name: OpenLab ("OL") and the "Mustang"
%   ChemStation build use 8-byte doubles; other GC builds use 4-byte floats.
%
% Intercept note (version 130):
%   chromConverter reads an 8-byte "intercept" double at offset 4110 for
%   version 130 and adds it. That is a bug: offsets 4110/4114 are two separate
%   int32 fields (an unidentified value, and the minimum intensity), so the
%   8-byte read yields a denormal (~1e-311) which is numerically a no-op.
%   Neither rainbow nor OpenChrom uses any intercept for this version, so this
%   importer uses intercept = 0 for version 130. For the other versions the
%   intercept is a genuine stored double and is applied.
%
% Input Arguments:
%   filePath - (1,1) string
%       Path to the .CH file.
%
% Name-Value Arguments:
%   ApplyScaling - (1,1) logical
%       When true (default), convert the raw stored signal to physical units
%       as intensity = raw * scalingFactor + intercept, in the unit reported
%       by signal.unit. For the validated diode-array file that scaling
%       factor happens to be exactly 1000/2^21 (a 21-bit ADC over a 1000 mAU
%       full scale) and the unit is "mAU"; other detectors store their own
%       factor, intercept and unit in the same header fields. When false, the
%       raw decoded values are returned and scalingFactor/intercept are left
%       for the caller to apply.
%
% Output Arguments:
%   fileContent - (1,1) struct with the fields below. Every field is always present;
%   a field the file's version does not provide, or provides but leaves empty,
%   is set to a missing/NaN value of its normal type, so that structs returned
%   by this function always share one schema and can be concatenated into an
%   array.
%
%       file - struct
%           .name    (1,1) string   file name
%           .path    (1,1) string   full path
%           .bytes   (1,1) double   file size in bytes
%           .type    (1,1) string   file-type string, e.g. "LC DATA FILE"
%           .version (1,1) string   version tag, e.g. "130"
%       sample - struct
%           .name         (1,1) string
%           .description  (1,1) string  only stored by versions 8/81; missing
%                                       otherwise
%           .vial         (1,1) double  autosampler vial number
%           .vialPosition (1,1) string  vial position as text, e.g. "D1B-A4";
%                                       only stored by version 130
%           .sequence     (1,1) double  sequence-line number
%           .replicate    (1,1) double  injection/replicate number
%       method - struct
%           .name        (1,1) string   method file, e.g. "OXI2B.M"
%           .operator    (1,1) string   operator/user name
%           .datetime    (1,1) datetime acquisition time; NaT if the raw text
%                                       could not be parsed
%           .rawDatetime (1,1) string   the raw header date text
%       instrument - struct
%           .name              (1,1) string  detector/inlet identifier
%           .inlet             (1,1) string  inlet/technique label; missing for
%                                            versions 179/181
%           .signalDescription (1,1) string  e.g. "DAD1A, Sig=254,16 Ref=..";
%                                            missing for versions 8/81
%           .firmwareRevision  (1,1) string  acquiring module firmware, e.g.
%                                            "B.07.10 [0004]"; only versions
%                                            30/130
%       software - struct
%           .name        (1,1) string  acquisition software, e.g.
%                                      "Asterix ChemStation"; missing for 8/81
%           .revision    (1,1) string  e.g. "Rev. C.01.10 [287]"; only 30/130
%       signal - struct (single chromatogram channel)
%           .retentionTime  (:,1) double   minutes since start of run
%           .intensity      (:,1) double   one value per retention time
%           .unit           (1,1) string   unit read from the file header,
%                                          e.g. "mAU" for a DAD/MWD channel
%           .minIntensity   (1,1) double   header-stored minimum; NaN unless
%                                          version 130
%           .maxIntensity   (1,1) double   header-stored maximum; NaN unless
%                                          version 130
%           .scalingFactor  (1,1) double   header scaling factor
%           .scalingApplied (1,1) logical  whether .intensity was scaled
%
%   For version 130, minIntensity/maxIntensity (header offsets 4114/4118) are
%   verified against the decoded signal on import; a mismatch raises a warning.
%
% This file is self-contained: the helper functions it uses (readAgilentString
% and parseDateTimeText) are included below as local functions.

arguments
    filePath (1,1) string {mustBeFile}
    options.ApplyScaling (1,1) logical = true
end

fileInfo = dir(filePath);

% Numeric header fields are big-endian; text is length-prefixed, either ASCII
% (legacy versions) or UTF-16LE (versions 130/179/181).
fid = fopen(filePath, "r", "b");
if fid < 0
    error("importAgilentCH:cannotOpen", "Could not open file: %s", filePath);
end
cleanup = onCleanup(@() fclose(fid));

% --- File version tag (offset 0: length-prefixed ASCII, e.g. "130") -------
versionTag = readAgilentString(fid, 0, Encoding="ascii");
layout = getLayout(versionTag, filePath);
enc = layout.encoding;

% Resolve the version-179 sub-form (4-byte float vs 8-byte double).
if layout.decoder == "doubleArray"
    layout.doubleArrayBytes = resolve179ByteWidth(fid);
end

% --- Metadata strings -----------------------------------------------------
fileType     = readOptionalString(fid, layout.fileType,     enc);
sampleName   = readOptionalString(fid, layout.sampleName,   enc);
description  = readOptionalString(fid, layout.description,  enc);
operator     = readOptionalString(fid, layout.operator,     enc);
methodName   = readOptionalString(fid, layout.method,       enc);
detector     = readOptionalString(fid, layout.detector,     enc);
inletName    = readOptionalString(fid, layout.inlet,        enc);
units        = readOptionalString(fid, layout.units,        enc);
signalDesc   = readOptionalString(fid, layout.signal,       enc);
software     = readOptionalString(fid, layout.software,     enc);
firmwareRevision = readOptionalString(fid, layout.firmwareRevision, enc);
softwareRevision = readOptionalString(fid, layout.softwareRevision, enc);
vialPosition = readOptionalString(fid, layout.vialPosition, enc);

rawDate = readOptionalString(fid, layout.date, enc);
methodDatetime = parseDateTimeText(rawDate);

% --- Sample sequence / vial / replicate (big-endian int16) ----------------
fseek(fid, 252, "bof");
sequenceLine = fread(fid, 1, "int16", 0, "b");
vial = fread(fid, 1, "int16", 0, "b");
replicate = fread(fid, 1, "int16", 0, "b");

% --- Data start -----------------------------------------------------------
if layout.decoder == "doubleArray"
    dataStart = 6144;   % fixed for the uncompressed-double form
else
    % Word pointer at offset 264 (byte = (value - 1) * 512).
    fseek(fid, 264, "bof");
    dataStart = (fread(fid, 1, "int32", 0, "b") - 1) * 512;
end

% --- Retention-time bounds (in ms; int32 or float32 depending on version) -
fseek(fid, layout.timeStart, "bof");
tStart = double(fread(fid, 1, layout.timeFormat, 0, "b")) / 60000;
fseek(fid, layout.timeEnd, "bof");
tEnd   = double(fread(fid, 1, layout.timeFormat, 0, "b")) / 60000;

% --- Calibration ----------------------------------------------------------
fseek(fid, layout.scalingFactor, "bof");
scalingFactor = fread(fid, 1, "float64", 0, "b");

% Version-8 scaling toggle: a small header code overrides the factor.
if ~isnan(layout.scalingToggle)
    fseek(fid, layout.scalingToggle, "bof");
    toggle = fread(fid, 1, "int32", 0, "b");
    if ismember(toggle, [1 2 3])
        scalingFactor = 1.33321110047553;
    end
end

intercept = 0;
if ~isnan(layout.intercept)
    fseek(fid, layout.intercept, "bof");
    value = fread(fid, 1, "float64", 0, "b");
    % Guard against a bogus/denormal value (see the intercept note above).
    if isfinite(value) && abs(value) > realmin("double")
        intercept = value;
    end
end

% --- Header-stored raw minimum / maximum (used as an integrity check) -----
rawMin = readOptionalInt32(fid, layout.minIntensity);
rawMax = readOptionalInt32(fid, layout.maxIntensity);

% --- Data body ------------------------------------------------------------
switch layout.decoder
    case "delta"
        intensity = decodeDelta(fid, dataStart, fileInfo.bytes);
    case "doubleDelta"
        intensity = decodeDoubleDelta(fid, dataStart, fileInfo.bytes);
    case "doubleArray"
        intensity = decodeDoubleArray(fid, dataStart, layout.doubleArrayBytes);
end
retentionTime = linspace(tStart, tEnd, numel(intensity))';

% Verify the decoded signal against the min/max recorded in the header.
if ~isnan(rawMin) && ~isnan(rawMax) && ~isempty(intensity)
    if min(intensity) ~= rawMin || max(intensity) ~= rawMax
        warning("importAgilentCH:minMaxMismatch", ...
            join(["Decoded signal range [%g %g] does not match the range" ...
                  "recorded in the header [%g %g] for %s. The decode may be" ...
                  "incorrect."]), ...
            min(intensity), max(intensity), rawMin, rawMax, filePath);
    end
end

scalingApplied = false;
minIntensity = rawMin;
maxIntensity = rawMax;
if options.ApplyScaling
    intensity = intensity * scalingFactor + intercept;
    minIntensity = rawMin * scalingFactor + intercept;
    maxIntensity = rawMax * scalingFactor + intercept;
    scalingApplied = true;
end

% --- Assemble output ------------------------------------------------------
fileContent = struct();
fileContent.file = struct("name", string(fileInfo.name), ...
                   "path", string(filePath), ...
                   "bytes", fileInfo.bytes, ...
                   "type", fileType, ...
                   "version", string(versionTag));
fileContent.sample = struct("name", sampleName, ...
                     "description", description, ...
                     "vial", vial, ...
                     "vialPosition", vialPosition, ...
                     "sequence", sequenceLine, ...
                     "replicate", replicate);
fileContent.method = struct("name", methodName, ...
                     "operator", operator, ...
                     "datetime", methodDatetime, ...
                     "rawDatetime", rawDate);
fileContent.instrument = struct("name", detector, ...
                         "inlet", inletName, ...
                         "signalDescription", signalDesc, ...
                         "firmwareRevision", firmwareRevision);
fileContent.software = struct("name", software, ...
                       "revision", softwareRevision);
fileContent.signal = struct("retentionTime", retentionTime, ...
                     "intensity", intensity, ...
                     "unit", units, ...
                     "minIntensity", minIntensity, ...
                     "maxIntensity", maxIntensity, ...
                     "scalingFactor", scalingFactor, ...
                     "scalingApplied", scalingApplied);

end


function layout = getLayout(versionTag, filePath)
%GETLAYOUT Byte offsets and decoding options for a given .CH version.
%   Offsets are zero-based (as used by fseek 'bof'); NaN means the version does
%   not provide that field. Offsets follow chromConverter's
%   get_agilent_offsets(); the version-130 min/max intensity offsets
%   (4114/4118) were additionally verified byte-for-byte against the datasets.
%
%   Output-field offsets ('detector'->instrument.name, 'inlet'->instrument.inlet,
%   'firmwareRevision', 'softwareRevision') are mapped so the returned struct is
%   consistent across versions even though Agilent reuses header slots
%   differently between the legacy and 130-family layouts.

switch versionTag
    case {"8", "81"}
        % Legacy FID/ADC. ASCII strings; has a sample description; no software
        % or signal-descriptor fields.
        layout = baseLayout( ...
            "encoding",         "ascii", ...
            "fileType",         4, ...
            "sampleName",       24, ...
            "description",      86, ...
            "operator",         148, ...
            "date",             178, ...
            "detector",         208, ...
            "inlet",            218, ...
            "method",           228, ...
            "units",            580, ...
            "intercept",        636, ...
            "scalingFactor",    644);
        if versionTag == "8"
            layout.decoder = "delta";
            layout.timeFormat = "int32";
            layout.scalingToggle = 542;   % version-8 only
        else
            layout.decoder = "doubleDelta";
            layout.timeFormat = "float32";
        end
    case "30"
        % Legacy DAD/MWD/ELSD. ASCII strings; software + signal descriptor.
        layout = baseLayout( ...
            "encoding",         "ascii", ...
            "fileType",         4, ...
            "sampleName",       24, ...
            "operator",         148, ...
            "date",             178, ...
            "detector",         208, ...
            "inlet",            218, ...
            "method",           228, ...
            "software",         322, ...
            "firmwareRevision", 355, ...   % unverified for v30 (see note below)
            "softwareRevision", 405, ...
            "units",            580, ...
            "signal",           596, ...
            "intercept",        636, ...
            "scalingFactor",    644, ...
            "decoder",          "delta", ...
            "timeFormat",       "int32");
    case "130"
        % Modern DAD/MWD/ELSD. UTF-16LE strings; software, firmware, signal,
        % vial position, and header min/max.
        layout = baseLayout( ...
            "encoding",         "utf16le", ...
            "fileType",         347, ...
            "sampleName",       858, ...
            "operator",         1880, ...
            "date",             2391, ...
            "detector",         2492, ...   % holds an inlet code ("GCI") here
            "inlet",            2533, ...   % technique ("LC")
            "method",           2574, ...
            "software",         3089, ...
            "firmwareRevision", 3601, ...   % module firmware "B.07.10 [0004]"
            "softwareRevision", 3802, ...
            "units",            4172, ...
            "signal",           4213, ...
            "vialPosition",     4055, ...   % same offset as .UV 131
            "minIntensity",     4114, ...
            "maxIntensity",     4118, ...
            "intercept",        NaN, ...    % see the intercept note in the help
            "scalingFactor",    4732, ...
            "decoder",          "delta", ...
            "timeFormat",       "int32");
    case {"179", "181"}
        % Modern FID/ADC. UTF-16LE strings; software name and signal, but no
        % firmware/description/vial-position/min-max.
        layout = baseLayout( ...
            "encoding",         "utf16le", ...
            "fileType",         347, ...
            "sampleName",       858, ...
            "operator",         1880, ...
            "date",             2391, ...
            "detector",         2492, ...
            "method",           2574, ...
            "software",         3089, ...
            "units",            4172, ...
            "signal",           4213, ...
            "intercept",        4724, ...
            "scalingFactor",    4732);
        if versionTag == "181"
            layout.decoder = "doubleDelta";
        else
            layout.decoder = "doubleArray";
        end
        layout.timeFormat = "float32";
    otherwise
        error("importAgilentCH:unsupportedVersion", ...
            join(["Unsupported .CH file version tag '%s' in %s. Supported" ...
                  "versions are 8, 30, 81, 130, 179 and 181."]), ...
            versionTag, filePath);
end

end


function layout = baseLayout(varargin)
%BASELAYOUT Build a layout struct with every field defaulted, then overridden.
%   Ensures all layouts expose the same field set (unset offsets = NaN) so the
%   importer body can read them uniformly.
defaults = struct( ...
    "encoding",         "utf16le", ...
    "fileType",         NaN, ...
    "sampleName",       NaN, ...
    "description",      NaN, ...
    "operator",         NaN, ...
    "date",             NaN, ...
    "detector",         NaN, ...
    "inlet",            NaN, ...
    "method",           NaN, ...
    "software",         NaN, ...
    "firmwareRevision", NaN, ...
    "softwareRevision", NaN, ...
    "units",            NaN, ...
    "signal",           NaN, ...
    "vialPosition",     NaN, ...
    "minIntensity",     NaN, ...
    "maxIntensity",     NaN, ...
    "intercept",        NaN, ...
    "scalingFactor",    NaN, ...
    "scalingToggle",    NaN, ...
    "timeStart",        282, ...
    "timeEnd",          286, ...
    "timeFormat",       "int32", ...
    "decoder",          "delta");
layout = defaults;
for k = 1:2:numel(varargin)
    layout.(varargin{k}) = varargin{k+1};
end
end


function bytesForm = resolve179ByteWidth(fid)
%RESOLVE179BYTEWIDTH Choose 4-byte vs 8-byte doubles for a version-179 file.
%   OpenLab ("OL") files and the "Mustang" ChemStation build use 8-byte
%   doubles; other GC builds use 4-byte floats.
typeCode = readAgilentString(fid, 347, Encoding="utf16le");
typeCode = extractBefore(typeCode + "  ", 3);
if typeCode == "OL"
    bytesForm = 8;
else
    software = readAgilentString(fid, 3089, Encoding="utf16le");
    firstWord = extractBefore(software + " ", " ");
    bytesForm = 4;
    if firstWord == "Mustang"
        bytesForm = 8;
    end
end
end


function value = readOptionalInt32(fid, offset)
%READOPTIONALINT32 Read a big-endian int32, or NaN when the offset is absent.
if isnan(offset)
    value = NaN;
    return
end
fseek(fid, offset, "bof");
value = fread(fid, 1, "int32", 0, "b");
if isempty(value)
    value = NaN;
end
end


function str = readOptionalString(fid, offset, encoding)
%READOPTIONALSTRING Read a string field, or missing when the offset is absent.
%   IDENTICAL COPY - this local function also exists in the other Agilent
%   importers (importAgilentUV/MS/CH). Keep the copies in sync.
if isnan(offset)
    str = string(missing);
else
    str = readAgilentString(fid, offset, Encoding=encoding);
end
end


function signal = decodeDelta(fid, dataStart, fileBytes)
%DECODEDELTA Decode the delta-compressed .CH signal (versions 8, 30, 130).
%   The body is a sequence of segments. Each segment starts with a marker byte
%   0x10 followed by a length byte giving the number of values in it. Values
%   are int16 deltas accumulated onto a running total; a delta of -32768
%   (0x8000) is an escape meaning the running total is instead reset from the
%   next 4-byte integer. Decoding stops at the first non-0x10 marker or EOF.

signal = zeros(ceil((fileBytes - dataStart) / 2), 1);
index = 1;
segmentCarry = 0;   % running total carried across segments

fseek(fid, dataStart, "bof");
while ftell(fid) < fileBytes
    marker = fread(fid, 1, "uint8", 0, "b");
    if isempty(marker) || marker ~= 16      % 0x10 starts a segment
        break
    end

    segmentLength = fread(fid, 1, "uint8", 0, "b");
    if isempty(segmentLength)
        break
    end

    accumulator = segmentCarry;
    for i = 1:segmentLength
        delta = fread(fid, 1, "int16", 0, "b");
        if isempty(delta)
            break
        end
        if delta ~= -32768
            accumulator = accumulator + delta;
        else
            accumulator = fread(fid, 1, "int32", 0, "b");
        end
        signal(index) = accumulator;
        index = index + 1;
    end
    segmentCarry = accumulator;
end

signal = signal(1:index-1);

end


function signal = decodeDoubleDelta(fid, dataStart, fileBytes)
%DECODEDOUBLEDELTA Decode a double-delta-compressed signal (versions 81, 181).
%   UNTESTED - ported from chromConverter's decode_double_delta (which credits
%   the Chromatography Toolbox). A running first and second difference are
%   accumulated from int16 values; the sentinel 32767 introduces a full
%   6-byte absolute value (int16 high word * 2^32 + uint32 low word).

signal = zeros(ceil((fileBytes - dataStart) / 2), 1);
count = 1;
b1 = 0;   % running value
b2 = 0;   % running first difference

fseek(fid, dataStart, "bof");
while ftell(fid) < fileBytes
    b3 = fread(fid, 1, "int16", 0, "b");
    if isempty(b3)
        break
    end
    if b3 ~= 32767
        b2 = b2 + b3;
        b1 = b1 + b2;
    else
        high = fread(fid, 1, "int16", 0, "b");
        low  = fread(fid, 1, "uint32", 0, "b");
        if isempty(high) || isempty(low)
            break
        end
        b1 = high * 4294967296 + low;
        b2 = 0;
    end
    signal(count) = b1;
    count = count + 1;
end

signal = signal(1:count-1);

end


function signal = decodeDoubleArray(fid, dataStart, byteWidth)
%DECODEDOUBLEARRAY Decode uncompressed little-endian reals (version 179).
%   UNTESTED - ported from chromConverter's decode_double_array_*. The 8-byte
%   form is a plain array of float64 values. The 4-byte form stores interleaved
%   pairs of float32 and the signal is every second value.

fseek(fid, dataStart, "bof");
if byteWidth == 8
    signal = fread(fid, Inf, "float64", 0, "l");
else
    raw = fread(fid, Inf, "float32", 0, "l");
    signal = raw(2:2:end);
end
signal = signal(:);

end


% ========================================================================
% Local helper functions
% ------------------------------------------------------------------------
% These two functions are duplicated verbatim in importAgilentUV,
% importAgilentMS and importAgilentCH so that each importer is a single,
% self-contained file. If you fix or extend one copy, update the other two.
% ========================================================================

function str = readAgilentString(fid, offset, options)
%READAGILENTSTRING Read a length-prefixed string from an Agilent data file.
%   IDENTICAL COPY - this local function also exists, byte-for-byte, in the
%   other Agilent importers (importAgilentUV/MS/CH). Keep the copies in sync.
%
%   Agilent files store text fields "Pascal style": a uint8 length followed by
%   that many characters, either one byte each (ASCII, older formats) or two
%   bytes each (UTF-16LE, newer formats). Returns "" for an empty/absent field.
arguments
    fid (1,1) double
    offset (1,1) double {mustBeNonnegative}
    options.Encoding (1,1) string {mustBeMember(options.Encoding, ["ascii", "utf16le"])} = "utf16le"
end
if fseek(fid, offset, "bof") ~= 0
    str = "";
    return
end
n = fread(fid, 1, "uint8");
if isempty(n) || n == 0
    str = "";
    return
end
switch options.Encoding
    case "ascii"
        raw = fread(fid, n, "uint8=>char", 0, "l")';
    case "utf16le"
        % Each character is stored as a little-endian uint16 code unit.
        raw = fread(fid, n, "uint16=>char", 0, "l")';
end
str = string(strtrim(deblank(raw)));
end


function dt = parseDateTimeText(dateTimeText)
%PARSEDATETIMETEXT Parse an Agilent header date/time string into a datetime.
%   IDENTICAL COPY - this local function also exists, byte-for-byte, in the
%   other Agilent importers (importAgilentUV/MS/CH). Keep the copies in sync.
%
%   Always returns a datetime: NaT on empty/missing input or when no known
%   format matches. Format is "preserveinput" and the input's time zone is
%   preserved (offset-bearing .MS dates come back zoned; .UV/.CH unzoned).
arguments
    dateTimeText (1,1) string = ""
end
dt = NaT;
if ismissing(dateTimeText)
    return
end
% Trim and collapse internal whitespace (some formats have double spaces).
txt = regexprep(strtrim(char(dateTimeText)), '\s+', ' ');
if isempty(txt)
    return
end
% Candidate formats, most specific first (LDML datetime patterns).
formats = { ...
    'dd MMM yy hh:mm a xxxx', ...   % 15 Jan 24 12:10 pm +0100
    'dd MMM yyyy hh:mm a xxxx', ... % 15 Jan 2024 12:10 pm +0100
    'dd MMM yy hh:mm a', ...        % 15 Jan 24 12:10 pm
    'dd-MMM-yy, HH:mm:ss', ...      % 15-Jan-24, 12:10:26
    'dd-MMM-yyyy, HH:mm:ss', ...    % 15-Jan-2024, 12:10:26
    'dd-MMM-yy HH:mm:ss', ...       % 15-Jan-24 12:10:26
    'MM/dd/yy hh:mm:ss a', ...      % 01/15/24 12:10:26 pm
    'MM/dd/yy HH:mm:ss', ...        % 01/15/24 12:10:26
    'MM/dd/yyyy HH:mm:ss'};         % 01/15/2024 12:10:26
for k = 1:numel(formats)
    fmt = formats{k};
    try
        if contains(fmt, 'x') || contains(fmt, 'X') || contains(fmt, 'Z')
            % Formats carrying a UTC offset need a TimeZone to parse.
            parsed = datetime(txt, InputFormat=fmt, Locale="en_US", ...
                TimeZone="local", Format="preserveinput");
        else
            parsed = datetime(txt, InputFormat=fmt, Locale="en_US", ...
                Format="preserveinput");
        end
    catch
        continue    % try the next format
    end
    if ~isnat(parsed)
        dt = parsed;
        return
    end
end
end
