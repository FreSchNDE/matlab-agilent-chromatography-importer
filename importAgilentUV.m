function fileContent = importAgilentUV(filePath, options)
%IMPORTAGILENTUV Import an Agilent ChemStation/OpenLab DAD .UV file.
%
% Syntax:
%   fileContent = importAgilentUV(filePath)
%   fileContent = importAgilentUV(filePath, ApplyScaling=false)
%
% Description:
%   Reads a 3D UV/Vis diode-array-detector spectrum file (e.g. DAD1.UV) from
%   an Agilent .D data folder and returns its metadata and the full
%   absorbance matrix (retention time x wavelength).
%
%   The .UV file stores, for every retention-time point, a complete
%   absorbance spectrum measured simultaneously across the diode array. The
%   individual DAD1A.ch / DAD1B.ch chromatograms are bandwidth-averaged,
%   reference-corrected single-wavelength views derived from this same cube;
%   this file is the raw source for all of them.
%
%   Supported file versions
%   -----------------------
%   Detected from the version tag at the start of the file (and, for "131",
%   the 2-character type code at offset 347):
%     "131" + "LC" : OpenLab CDS "LC" variant, delta-encoded    [VALIDATED]
%     "131" + "OL" : OpenLab CDS "OL" variant, 8-byte doubles   [UNTESTED]
%     "31"         : legacy ChemStation, delta-encoded          [UNTESTED]
%   Only 131_LC files were available to test; the "131_OL" (uncompressed
%   doubles) and legacy "31" paths are implemented from the chromConverter /
%   rainbow format references but are marked UNTESTED. An unrecognized version
%   tag is a hard error.
%
%   The binary layout (header offsets, wavelength encoding, delta scheme and
%   scaling factor) was validated byte-for-byte against the parser in the
%   chromConverter R package (read_chemstation_uv, adapted from the rainbow
%   project), which in turn derives from the Chromatography Toolbox
%   (c) James Dillon 2014. See tests/TestImportAgilent.m for the automated
%   comparison.
%
% Input Arguments:
%   filePath - (1,1) string
%       Path to the .UV file.
%
% Name-Value Arguments:
%   ApplyScaling - (1,1) logical
%       When true (default), multiply the raw integer counts by the header
%       scaling factor to obtain the physical unit reported in signal.unit
%       (typically "mAU" for a UV/Vis DAD). When false, the raw stored counts
%       are returned unscaled. The scaling factor is returned in either case,
%       so the operation is reversible.
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
%           .version (1,1) string   detected version key, e.g. "131_LC"
%       sample - struct
%           .name         (1,1) string
%           .vial         (1,1) double  autosampler vial number
%           .vialPosition (1,1) string  vial position as text, e.g. "D1B-A4";
%                                       missing for legacy version 31
%           .sequence     (1,1) double  sequence-line number
%           .replicate    (1,1) double  injection/replicate number
%       method - struct
%           .name        (1,1) string   method file, e.g. "OXI2B.M"
%           .operator    (1,1) string   operator/user name
%           .datetime    (1,1) datetime acquisition time; NaT if the raw text
%                                       could not be parsed
%           .rawDatetime (1,1) string   the raw header date text
%       instrument - struct
%           .name              (1,1) string  detector, e.g. "DAD1"
%           .inlet             (1,1) string  technique/inlet, e.g. "LC"
%           .signalDescription (1,1) string  e.g. "DAD1I, DAD: Spectrum";
%                                            missing for legacy version 31
%       signal - struct (UV/Vis DAD spectrum)
%           .retentionTime  (:,1) double   minutes since start of run
%           .wavelength     (1,:) double   wavelengths in nm
%           .absorbance     (:,:) double   (retentionTime x wavelength)
%           .unit           (1,1) string   unit read from the file header,
%                                          e.g. "mAU"
%           .scalingFactor  (1,1) double   header scaling factor
%           .scalingApplied (1,1) logical  whether .absorbance was scaled
%
% This file is self-contained: the helper functions it uses (readAgilentString
% and parseDateTimeText) are included below as local functions.

arguments
    filePath (1,1) string {mustBeFile}
    options.ApplyScaling (1,1) logical = true
end

fileInfo = dir(filePath);

% The file mixes endianness: most of the header numerics and the entire data
% body are little-endian, while a few header fields (number of scans, scaling
% factor) are big-endian. Open little-endian by default and override per read.
fid = fopen(filePath, "r", "l");
if fid < 0
    error("importAgilentUV:cannotOpen", "Could not open file: %s", filePath);
end
cleanup = onCleanup(@() fclose(fid));

% --- File version tag (offset 0: length-prefixed ASCII, e.g. "131" / "31") ---
versionTag = readAgilentString(fid, 0, Encoding="ascii");

% Version "131" additionally carries a length-prefixed type string at offset
% 347 (e.g. "LC Diode Array..."); its first two characters ("LC" or "OL")
% select the data encoding.
switch versionTag
    case "131"
        typeCode = readAgilentString(fid, 347, Encoding="utf16le");
        typeCode = extractBefore(typeCode + "  ", 3);  % first two characters
        versionKey = "131_" + upper(typeCode);
    case "31"
        versionKey = "31";
    otherwise
        error("importAgilentUV:unsupportedVersion", ...
            "Unsupported .UV file version tag '%s' in %s.", versionTag, filePath);
end

layout = getLayout(versionKey, filePath);

% --- Metadata strings -----------------------------------------------------
fileType     = readOptionalString(fid, layout.fileType,   layout.encoding);
sampleName   = readAgilentString(fid, layout.sampleName, Encoding=layout.encoding);
operator     = readAgilentString(fid, layout.operator,   Encoding=layout.encoding);
methodName   = readAgilentString(fid, layout.method,     Encoding=layout.encoding);
detectorName = readAgilentString(fid, layout.detector,   Encoding=layout.encoding);
inletName    = readOptionalString(fid, layout.instrument, layout.encoding);
signalDesc   = readOptionalString(fid, layout.signal,     layout.encoding);
units        = readAgilentString(fid, layout.units,      Encoding=layout.encoding);

% Autosampler vial position as text (e.g. a well-plate coordinate like
% "D1B-A4"). Only stored by the 131 family; missing for legacy version 31.
% This complements the numeric 'sample.vial' read below.
vialPosition = readOptionalString(fid, layout.vialPosition, layout.encoding);

rawDate = readAgilentString(fid, layout.date, Encoding=layout.encoding);
methodDatetime = parseDateTimeText(rawDate);

% --- Header numerics ------------------------------------------------------
% Sample sequence line, vial number and replicate/injection number are stored
% as consecutive big-endian int16 at fixed offsets (shared with .MS/.CH).
fseek(fid, 252, "bof");
sequenceLine = fread(fid, 1, "int16", 0, "b");
vial = fread(fid, 1, "int16", 0, "b");
replicate = fread(fid, 1, "int16", 0, "b");

% Number of retention-time points (4-byte big-endian integer).
fseek(fid, layout.numTimes, "bof");
nTimes = fread(fid, 1, "int32", 0, "b");

% Scaling factor (8-byte big-endian double).
fseek(fid, layout.scalingFactor, "bof");
scalingFactor = fread(fid, 1, "float64", 0, "b");

% Wavelength axis: three little-endian int16 at dataStart+8, in units of
% 1/20 nm (start, end, step).
fseek(fid, layout.dataStart + 8, "bof");
waveInfo = fread(fid, 3, "int16", 0, "l");
lambdaStart = double(waveInfo(1)) / 20;
lambdaEnd   = double(waveInfo(2)) / 20;
lambdaStep  = double(waveInfo(3)) / 20;
wavelength = lambdaStart:lambdaStep:lambdaEnd;
nWave = numel(wavelength);

% --- Data body ------------------------------------------------------------
switch layout.decoder
    case "delta"
        [retentionTime, absorbance] = decodeDelta(fid, layout.dataStart, nTimes, nWave);
    case "array"
        [retentionTime, absorbance] = decodeArray(fid, layout.dataStart, nTimes, nWave);
end

scalingApplied = false;
if options.ApplyScaling
    absorbance = absorbance * scalingFactor;
    scalingApplied = true;
end

% --- Assemble output ------------------------------------------------------
fileContent = struct();
fileContent.file = struct("name", string(fileInfo.name), ...
                   "path", string(filePath), ...
                   "bytes", fileInfo.bytes, ...
                   "type", fileType, ...
                   "version", versionKey);
fileContent.sample = struct("name", sampleName, ...
                     "vial", vial, ...
                     "vialPosition", vialPosition, ...
                     "sequence", sequenceLine, ...
                     "replicate", replicate);
fileContent.method = struct("name", methodName, ...
                     "operator", operator, ...
                     "datetime", methodDatetime, ...
                     "rawDatetime", rawDate);
fileContent.instrument = struct("name", detectorName, ...
                         "inlet", inletName, ...
                         "signalDescription", signalDesc);
fileContent.signal = struct("retentionTime", retentionTime, ...
                     "wavelength", wavelength, ...
                     "absorbance", absorbance, ...
                     "unit", units, ...
                     "scalingFactor", scalingFactor, ...
                     "scalingApplied", scalingApplied);

end


function str = readOptionalString(fid, offset, encoding)
%READOPTIONALSTRING Read a string field whose offset may be NaN (absent).
%   IDENTICAL COPY - this local function also exists in the other Agilent
%   importers (importAgilentUV/MS/CH). Keep the copies in sync.
%
%   Returns a missing string (not "") when the field is absent, so callers can
%   distinguish "not in this format" from a genuinely empty stored string.
if isnan(offset)
    str = string(missing);
else
    str = readAgilentString(fid, offset, Encoding=encoding);
end
end


function layout = getLayout(versionKey, filePath)
%GETLAYOUT Byte offsets and decoding options for a given .UV version.
%   Offsets are zero-based (as used by fseek 'bof'). String offsets that the
%   version does not provide are NaN and yield "". Offsets follow
%   chromConverter's get_agilent_offsets(), extended with fields ('fileType',
%   'instrument', 'signal') verified byte-for-byte against these datasets.

switch versionKey
    case {"131_LC", "131_OL"}
        layout = struct( ...
            "encoding",      "utf16le", ...
            "fileType",      347, ...
            "sampleName",    858, ...
            "operator",      1880, ...
            "date",          2391, ...
            "detector",      2492, ...
            "instrument",    2533, ...
            "method",        2574, ...
            "units",         3093, ...
            "signal",        3136, ...
            "vialPosition",  4055, ...
            "numTimes",      278, ...
            "scalingFactor", 3085, ...
            "dataStart",     4096, ...
            "decoder",       "delta");
        if versionKey == "131_OL"
            layout.decoder = "array";
        end
    case "31"
        layout = struct( ...
            "encoding",      "ascii", ...
            "fileType",      4, ...
            "sampleName",    24, ...
            "operator",      148, ...
            "date",          178, ...
            "detector",      208, ...
            "instrument",    218, ...
            "method",        228, ...
            "units",         326, ...
            "signal",        NaN, ...
            "vialPosition",  NaN, ...
            "numTimes",      278, ...
            "scalingFactor", 318, ...
            "dataStart",     512, ...
            "decoder",       "delta");
    otherwise
        error("importAgilentUV:unsupportedVersion", ...
            "Unsupported .UV file version '%s' in %s.", versionKey, filePath);
end

end


function [retentionTime, signal] = decodeDelta(fid, dataStart, nTimes, nWave)
%DECODEDELTA Decode a delta-compressed .UV data body (versions 31, 131_LC).
%   Each retention-time record is a 22-byte header (4 discarded bytes, a
%   4-byte little-endian time in ms, 14 discarded bytes) followed by nWave
%   little-endian int16 deltas. A delta of -32768 (0x8000) is an escape: the
%   accumulator is instead reset from the following 4-byte little-endian
%   integer (absolute value).

timeMs = zeros(nTimes, 1);
signal = zeros(nTimes, nWave);

fseek(fid, dataStart, "bof");
for i = 1:nTimes
    fread(fid, 1, "int32", 0, "l");             % discard record marker
    timeMs(i) = fread(fid, 1, "int32", 0, "l"); % retention time in ms
    fseek(fid, 14, "cof");                       % discard 14 bytes

    accum = 0;
    for j = 1:nWave
        delta = fread(fid, 1, "int16", 0, "l");
        if delta == -32768
            accum = fread(fid, 1, "int32", 0, "l");
        else
            accum = accum + delta;
        end
        signal(i, j) = accum;
    end
end

retentionTime = timeMs / 60000;

end


function [retentionTime, signal] = decodeArray(fid, dataStart, nTimes, nWave)
%DECODEARRAY Decode an uncompressed .UV data body (version 131_OL).
%   Same 22-byte per-record header as the delta format, but the spectrum is
%   stored directly as nWave little-endian 8-byte doubles.

timeMs = zeros(nTimes, 1);
signal = zeros(nTimes, nWave);

fseek(fid, dataStart, "bof");
for i = 1:nTimes
    fread(fid, 1, "int32", 0, "l");             % discard record marker
    timeMs(i) = fread(fid, 1, "int32", 0, "l"); % retention time in ms
    fseek(fid, 14, "cof");                       % discard 14 bytes
    signal(i, :) = fread(fid, nWave, "float64", 0, "l")';
end

retentionTime = timeMs / 60000;

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
