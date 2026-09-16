function fileContent = importAgilentMS(filePath, options)
%IMPORTAGILENTMS Import an Agilent ChemStation/OpenLab MSD .MS file.
%
% Syntax:
%   fileContent = importAgilentMS(filePath)
%   fileContent = importAgilentMS(filePath, Precision=2)
%
% Description:
%   Reads a single mass-spectrometer data file (e.g. MSD1.MS) from an Agilent
%   .D data folder and returns its metadata together with the total ion
%   current (TIC) and the full ion-abundance matrix (retention time x m/z).
%
%   The binary layout was reverse-engineered with reference to the rainbow
%   project's format documentation and the chromConverter R package, which in
%   turn derive from the Chromatography Toolbox (c) James Dillon 2014, and
%   validated byte-for-byte against chromConverter's own output (see
%   tests/TestImportAgilent.m).
%
%   Supported file versions
%   -----------------------
%   Header version tag (first bytes of the file): "2" or "20". Two variants
%   exist, detected from the file-type string at offset 4:
%       "MSD Spectral File"  -> LC-MS (e.g. MSD1.MS)   [VALIDATED]
%       "GC / MS Data File"  -> GC-MS (e.g. DATA.MS)   [UNTESTED]
%   They share the same data body but differ in their header: the scan count
%   sits at a different offset AND endianness (LC: big-endian at 0x118; GC:
%   little-endian at 0x142), and their metadata strings live at different
%   offsets. Reading an LC file with the GC offsets yields garbage, so the
%   variant is detected rather than assumed; an unrecognised type string is a
%   hard error. GC-MS support is implemented from documentation only, since no
%   GC-MS files were available to test it; the sparse GC metadata reflects
%   undocumented header fields, but the signal data should decode correctly.
%
%   The data body itself is read by walking the scan segments sequentially
%   from the header-length pointer at 0x10A (documented as applying to both
%   variants), rather than the LC-only trailing scan directory; both
%   approaches were verified to give identical results on the LC files here.
%
% Input Arguments:
%   filePath - (1,1) string
%       Path to the .MS file.
%
% Name-Value Arguments:
%   Precision - (1,1) double
%       Number of decimal places the m/z axis is rounded to (default 3).
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
%           .type    (1,1) string   file-type string, e.g. "MSD Spectral File"
%           .version (1,1) string   header version tag, e.g. "2"
%       sample - struct
%           .name        (1,1) string
%           .description (1,1) string  "" if present but empty; missing if the
%                                      version does not store it (e.g. GC-MS)
%           .vial        (1,1) double  autosampler vial number (NaN for GC-MS)
%           .sequence    (1,1) double  sequence-line number (NaN for GC-MS)
%           .replicate   (1,1) double  injection/replicate number (NaN for GC-MS)
%       method - struct
%           .name        (1,1) string   acquisition method file, e.g. "OXI2B.M"
%           .operator    (1,1) string   operator/user name
%           .datetime    (1,1) datetime acquisition time; NaT if the raw text
%                                       could not be parsed or is absent
%           .rawDatetime (1,1) string   the raw header date text
%       instrument - struct
%           .name              (1,1) string  instrument/detector name
%           .inlet             (1,1) string  inlet label (missing for GC-MS)
%           .signalDescription (1,1) string  e.g. "MSD1, Initial Scan Range=..."
%       signal - struct (mass spectra over time)
%           .retentionTime (:,1) double  minutes since start of run
%           .mz            (1,:) double  mass-to-charge axis
%           .tic           (:,1) double  total ion current (one per scan)
%           .xic           (:,:) double  ion abundance (retentionTime x mz)
%
% This file is self-contained: the helper functions it uses (readAgilentString
% and parseDateTimeText) are included below as local functions.

arguments
    filePath (1,1) string {mustBeFile}
    options.Precision (1,1) double {mustBeNonnegative, mustBeInteger} = 3
end

fileInfo = dir(filePath);

% Most .MS numeric fields are big-endian (the GC scan count is the exception,
% handled explicitly below). Open big-endian.
fid = fopen(filePath, "r", "b");
if fid < 0
    error("importAgilentMS:cannotOpen", "Could not open file: %s", filePath);
end
cleanup = onCleanup(@() fclose(fid));

% --- File version tag (offset 0: length-prefixed ASCII, e.g. "2" / "20") ---
versionTag = readAgilentString(fid, 0, Encoding="ascii");
if ~ismember(versionTag, ["2", "20"])
    error("importAgilentMS:unsupportedVersion", ...
        "Unsupported .MS file version tag '%s' in %s.", versionTag, filePath);
end

% --- Detect the LC-MS / GC-MS variant from the file-type string -----------
fileType = readAgilentString(fid, 4, Encoding="ascii");
layout = getLayout(fileType, filePath);

% --- Metadata strings -----------------------------------------------------
enc = layout.encoding;
sampleName  = readOptionalString(fid, layout.sampleName,  enc);
description = readOptionalString(fid, layout.description, enc);
operator    = readOptionalString(fid, layout.operator,    enc);
methodName  = readOptionalString(fid, layout.method,      enc);
instrument  = readOptionalString(fid, layout.instrument,  enc);
inlet       = readOptionalString(fid, layout.inlet,       enc);
signalInfo  = readOptionalString(fid, layout.signal,      enc);

rawDate = readOptionalString(fid, layout.date, enc);
methodDatetime = parseDateTimeText(rawDate);

% --- Sample sequence / vial / replicate (big-endian int16) ----------------
% Only defined for the LC header layout; NaN for GC (undocumented there).
[sequenceLine, vial, replicate] = deal(NaN);
if ~isnan(layout.sequence)
    fseek(fid, layout.sequence, "bof");
    sequenceLine = fread(fid, 1, "int16", 0, "b");
    vial = fread(fid, 1, "int16", 0, "b");
    replicate = fread(fid, 1, "int16", 0, "b");
end

% --- Number of scans: different offset AND endianness per variant ---------
fseek(fid, layout.numTimes, "bof");
nScans = fread(fid, 1, "uint16", 0, layout.numTimesEndian);

% --- Data start: header length in shorts at 0x10A (byte = value*2 - 2) ----
% Documented by rainbow as applying to both LC and GC files.
fseek(fid, 266, "bof");
dataStart = fread(fid, 1, "uint16", 0, "b") * 2 - 2;

% --- Walk the scan segments ----------------------------------------------
[retentionTime, totalIntensity, mz, xic] = ...
    readSpectra(fid, dataStart, nScans, fileInfo.bytes, options.Precision, filePath);

% --- Assemble output ------------------------------------------------------
fileContent = struct();
fileContent.file = struct("name", string(fileInfo.name), ...
                   "path", string(filePath), ...
                   "bytes", fileInfo.bytes, ...
                   "type", fileType, ...
                   "version", versionTag);
fileContent.sample = struct("name", sampleName, ...
                     "description", description, ...
                     "vial", vial, ...
                     "sequence", sequenceLine, ...
                     "replicate", replicate);
fileContent.method = struct("name", methodName, ...
                     "operator", operator, ...
                     "datetime", methodDatetime, ...
                     "rawDatetime", rawDate);
fileContent.instrument = struct("name", instrument, ...
                         "inlet", inlet, ...
                         "signalDescription", signalInfo);
fileContent.signal = struct("retentionTime", retentionTime, ...
                     "mz", mz, ...
                     "tic", totalIntensity, ...
                     "xic", xic);

end


function layout = getLayout(fileType, filePath)
%GETLAYOUT Header offsets for the LC-MS and GC-MS variants of the .MS format.
%   Offsets are zero-based; NaN means the variant does not provide that field.
%   LC offsets are validated against real files. GC offsets come from the
%   rainbow project's format documentation and are UNVALIDATED.

if contains(fileType, "MSD Spectral", IgnoreCase=true)
    % LC-MS: length-prefixed ASCII metadata, scan count big-endian at 0x118.
    layout = struct( ...
        "variant",         "LC", ...
        "encoding",        "ascii", ...
        "sampleName",      24, ...
        "description",     86, ...
        "operator",        148, ...
        "date",            178, ...
        "instrument",      208, ...
        "inlet",           218, ...
        "method",          228, ...
        "signal",          320, ...   % "MSD1, Initial Scan Range=200.0-2000.0"
        "sequence",        252, ...   % then vial (254) and replicate (256)
        "numTimes",        280, ...   % 0x118
        "numTimesEndian",  "b");
elseif contains(fileType, "GC / MS", IgnoreCase=true) || contains(fileType, "GC/MS", IgnoreCase=true)
    % GC-MS: UNVALIDATED. Metadata is null-separated (UTF-16LE-style) and most
    % header fields are undocumented, so only what rainbow identifies is read.
    % The scan count is LITTLE-endian at 0x142 - reading it big-endian, or at
    % the LC offset, yields garbage.
    layout = struct( ...
        "variant",         "GC", ...
        "encoding",        "utf16le", ...
        "sampleName",      NaN, ...
        "description",     NaN, ...
        "operator",        NaN, ...
        "date",            NaN, ...
        "instrument",      448, ...   % 0x1C0, e.g. "5977B GCM"
        "inlet",           NaN, ...
        "method",          1126, ...  % 0x466, e.g. "Rt-bDEX-SE_mcminn.M"
        "signal",          NaN, ...
        "sequence",        NaN, ...
        "numTimes",        322, ...   % 0x142
        "numTimesEndian",  "l");
else
    error("importAgilentMS:unknownVariant", ...
        join(["Unrecognised .MS file-type string '%s' in %s. Expected" ...
              "'MSD Spectral File' (LC-MS) or 'GC / MS Data File' (GC-MS)." ...
              "Refusing to guess, because LC and GC store the scan count at" ...
              "different offsets and endianness."]), fileType, filePath);
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


function [retentionTime, ticSignal, mzAxis, xic] = readSpectra(fid, dataStart, nScans, fileBytes, precision, filePath)
%READSPECTRA Walk the scan segments and build a dense (scan x m/z) matrix.
%   Segments are stored back to back from dataStart, so they are read
%   sequentially rather than via the (LC-only) trailing scan directory. This
%   is what rainbow does and works for both the LC and GC variants.
%
%   Each segment is:
%       int16  length of the whole segment, in 2-byte words
%       int32  retention time in ms
%       ...    (18-byte header in total)
%       n x    [uint16 mass, uint16 abundance] pairs
%       ...    10-byte footer whose last int32 is the TIC
%
%   Masses are uint16 fixed-point scaled by 20 (m/z * 20, i.e. 0.05 m/z
%   steps); abundances are a 14-bit mantissa with a 2-bit base-8 exponent.

retentionTime = zeros(nScans, 1);
ticSignal = zeros(nScans, 1);
mzAll = [];
abundanceAll = [];
pairsPerScan = zeros(nScans, 1);

position = dataStart;
for i = 1:nScans
    if position + 4 > fileBytes
        error("importAgilentMS:truncated", ...
            "Ran past the end of %s at scan %d of %d.", filePath, i, nScans);
    end

    fseek(fid, position, "bof");
    segmentWords = fread(fid, 1, "int16", 0, "b");
    segmentBytes = segmentWords * 2;
    nPairs = (segmentWords - 18) / 2 + 2;
    pairsPerScan(i) = nPairs;

    retentionTime(i) = fread(fid, 1, "int32", 0, "b") / 60000;   % minutes

    % Mass values (uint16, with the 2-byte abundance interleaved).
    fseek(fid, position + 18, "bof");
    mzAll(end+1:end+nPairs) = fread(fid, nPairs, "uint16", 2, "b");

    % Abundance values (uint16, interleaved with the masses).
    fseek(fid, position + 20, "bof");
    abundanceAll(end+1:end+nPairs) = fread(fid, nPairs, "uint16", 2, "b");

    % TIC: the last int32 of the 10-byte segment footer.
    fseek(fid, position + segmentBytes - 4, "bof");
    ticSignal(i) = fread(fid, 1, "int32", 0, "b");

    position = position + segmentBytes;
end

% Decode abundance: 14-bit mantissa, 2-bit base-8 exponent.
abundanceAll = bitand(abundanceAll, 16383, "uint16") .* ...
    (8 .^ abs(bitshift(abundanceAll, -14, "uint16")));

% Decode mass: stored as a fixed-point integer scaled by 20 (i.e. m/z * 20),
% giving 0.05 m/z resolution. Divide by 20, then round to requested precision.
mzAll = mzAll ./ 20;
mzAll = round(mzAll .* 10^precision) ./ 10^precision;

mzAxis = unique(mzAll, "sorted");

% Map each scan's pairs into the dense matrix.
stop = cumsum(pairsPerScan);
start = [1; stop(1:end-1) + 1];

xic = zeros(nScans, numel(mzAxis));
[~, cols] = ismember(mzAll, mzAxis);
for i = 1:nScans
    idx = start(i):stop(i);
    xic(i, cols(idx)) = abundanceAll(idx);
end

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
