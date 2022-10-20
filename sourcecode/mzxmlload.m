function out = mzxmlload(filename, nvpArgs)
% MZXMLREAD Read an mzXML file into MATLAB as a structure.
%
%   OUT = MZXMLREAD(FILENAME) reads an mzXML file into MATLAB as a
%   structure. FILENAME is a character vector or string scalar containing a
%   file name, or a path and a file name, of a mzXML file that conforms to 
%   the mzXML 2.1 specification or earlier specifications. OUT is a MATLAB 
%   structure with the following fields: scan, index, mzXML. The fields scan
%   and index are placed into the first level of the output structure for 
%   improved access to the spectra data, the remainder of the mzXML document
%   tree is parsed in following the schema specifications. The mzXML 2.1 
%   specification can be found at:
%
%   http://sashimi.sourceforge.net/schema_revision/mzXML_2.1/Doc/mzXML_2.1_tutorial.pdf
%
%   MZXMLREAD(...,'LEVELS',LEVEL) specifies which msLevel of scans to
%   extract from the mzXML file.  LEVEL must be a positive integer or
%   vector of integers. Default is to extract all scan levels. If you
%   specify LEVELS, you cannot use SCANINDICES.
%
%   MZXMLREAD(...,'TIMERANGE', RANGE) specifies which range of time to
%   extract scans from the mzXML file.  RANGE is a two-element numeric
%   array, [START END].  START is a scalar that must fall between the
%   startTime and endTime attributes of msRun element.  END is a scalar
%   that must fall between START and the endTime attribute of msRun.  START
%   must be less then END.  Default is to extract all scans. If you specify
%   TIMERANGE, you cannot use SCANINDICES.
%
%   MZXMLREAD(...,'SCANINDICES',INDICES) specifies an index or a vector of
%   indices extract scans. INDICES must contain unique positive integers
%   less or equal that the number of scans in the file. For example, to
%   indicate a range of indices, use [START_IND:END_IND]. Default is to
%   extract all scans. If you specify SCANINDICES, you cannot use LEVELS or
%   TIMERANGE.
%
%   MZXMLREAD(...,'VERBOSE',T/F) show progress of reading the mzXML file.
%   Default is false.
%
%   Example:
%
%       out = mzxmlread('results.mzxml');
%       % view the first scan in the mzXML file:
%       m = out.scan(1).peaks.mz(1:2:end);
%       z = out.scan(1).peaks.mz(2:2:end);
%       stem(m,z,'marker','none')
%
%   Note that the file results.mzxml is not provided. Sample files can be
%   found at http://sashimi.sourceforge.net/repository.html.
%
%   If you receive any errors related to memory or java heap space, try
%   increasing your java heap space as described here:
%
%       https://www.mathworks.com/support/solutions/data/1-18I2C.html
%
%   See also JCAMPREAD, MZCDF2PEAKS, MZCDFINFO, MZCDFREAD, MZXML2PEAKS,
%   MZXMLINFO, TGSPCINFO, TGSPCREAD.

% Copyright 2006-2021 The MathWorks, Inc.

% EDIT: Some MS files contain scan indexing errors where the scans do not start at 1 and increase at a non-constant value "Looking at You Agilent". A fix is implemented at line 145: Replace Original Scan Numbers to start at 1 and increment by 1

arguments
    filename(1,1) string
    nvpArgs.Levels(:,1) double {mustBeNonnegative, mustBeInteger}
    nvpArgs.TimeRange(1,2) double
    nvpArgs.ScanIndices(:,1) double {mustBeNonnegative, mustBeInteger}
    nvpArgs.Verbose(1,1) logical = false;
end

% Conflicting NVPs. ScanIndices is not compatible with TimeRange and Levels
if isfield(nvpArgs, "ScanIndices") && any(isfield(nvpArgs, ["Levels", "TimeRange"]))
    error(message('bioinfo:mzxmlread:ConflictingNVPs'));
end

% Check if file exists
if ~exist(filename,'file')
     error(message('bioinfo:mzxmlread:FileNotFound', filename));
end

% Grab first few lines of file to check if XML file and mzXML format
fid = fopen(filename,'rt');
str = fread(fid,120,'*char')';
fclose(fid);

% Check if XML file
if isempty(regexp(str,'<\?xml','once'))
    error(message('bioinfo:mzxmlread:missingXMLdeclaration', filename, filename));
end

% Check if mzXML format
if isempty(regexp(str,'<mzXML|<msRun','once'))
    error(message('bioinfo:mzxmlread:notValidMZXMLFile', filename));
end

% Check if filename contains full path, if not get it
if isempty(regexp(filename,filesep,'once'))
    filename = which(filename);
end

% Check if "ScanIndices" NVP exceeds number of scans
mzinfo = mzxmlinfo(filename);
if isfield(nvpArgs, "ScanIndices")
    if any(nvpArgs.ScanIndices > mzinfo.NumberOfScans)
        error(message('bioinfo:mzxmlread:InvalidScanIndices', 1, mzinfo.NumberOfScans));
    end
    scanIndices = nvpArgs.ScanIndices;
else
    scanIndices = 1:mzinfo.NumberOfScans;
end

% Read mzxml data.
if nvpArgs.Verbose
    disp("Starting to parse document...")
end
out = readstruct(filename, 'FileType', 'xml', "AttributeSuffix", "");

% Rename "Text" Field to "value" for backwards compatibility
if isfield(out, "index") && isfield(out.index, "offset") && isfield(out.index.offset, "Text")
    out.index.offset = renameStructField(out.index.offset, 'Text', 'value');
end

% Some mzxml files are improperly formatted without a "mzXML" top-level
% element, and instead have "msRun" as the top-level. Reformat so that the
% top level only contains mzXML fields.
outFields = string(fieldnames(out));
mzxmlFields = ["xmlns", "xmlns_xsi", "xsi_schemaLocation", "msRun", "index", "indexOffset", "sha1"];
msRunFieldLoc = find(~ismember(outFields, mzxmlFields))';
for iF = msRunFieldLoc
    out.msRun.(outFields(iF)) = out.(outFields(iF));
    out = rmfield(out, outFields(iF));
end

% Reorganize Struct hierarchy to fit the mzXML 2.1 specification found at:
% http://sashimi.sourceforge.net/schema_revision/mzXML_2.1/Doc/mzXML_2.1_tutorial.pdf
mzxml21Fields = ["xmlns", "xmlns_xsi", "xsi_schemaLocation", "msRun", "indexOffset", "sha1"];
mzxmlFieldsLoc = find(isfield(out, mzxml21Fields));
for iF = mzxmlFieldsLoc
    out.mzXML.(mzxml21Fields(iF)) = out.(mzxml21Fields(iF));
    out = rmfield(out, mzxml21Fields(iF));
    if strcmp(mzxml21Fields(iF), "msRun") && isfield(out.mzXML.msRun, "scan")
        out.scan = out.mzXML.msRun.scan;
        out.mzXML.msRun = rmfield(out.mzXML.msRun, "scan");
    end
end
% Replace Original Scan Numbers to start at 1 and increment by 1
scanNum = 1:mzinfo.NumberOfScans;
scanNum = num2cell(scanNum);
[out.scan(1:end).num] = scanNum{:};

% Format Scan Field. Move all scan fields to top level in the struct
% hierarchy and format so that all scan structs have the same fields
if isfield(out, "scan")
    if nvpArgs.Verbose
        disp("Formatting scan substructure...")
    end
    
    % Determine top-level fields for scans, including scans within scans
    topLevelFields = getTopLevelFieldNamess(out.scan);
    topLevelFields = union(fieldnames(formatScan()), topLevelFields);

    % Parse all scans and store in cell array. Add extra topLevelFields.
    % Concatenate formatted scans
    formattedScans = cell(numel(out.scan), 1);
    for iS = 1:numel(out.scan)
        formattedScans{iS} = parseScan(out.scan(iS), scanIndices);
        if ~isempty(formattedScans{iS})
            formattedScans{iS} = addExtraFields(formattedScans{iS}, topLevelFields);
        end
    end
    formattedScans = vertcat(formattedScans{:});

    % write to out.scan
    if isempty(formattedScans)
        out.scan = formatScan();
    else
        out.scan = formattedScans;
        [~, numOrder] = sort([out.scan.num]);
        out.scan = out.scan(numOrder);
    end
    
    % The inner text of an xml element is given the fieldname "Text" by
    % readstruct. Decode scan.peaks.Text and rename to scan.peaks.mz
    for iScan = 1:numel(out.scan)
        out.scan(iScan).peaks.mz = processPeaks(...
            typecast(matlab.net.base64decode(out.scan(iScan).peaks.Text), 'int8'),...
            out.scan(iScan).peaks.precision)';
        out.scan(iScan).peaks = rmfield(out.scan(iScan).peaks, "Text");

        % Rename precursorMz.Text -> precursorMz.value
        if isfield(out.scan(iScan).precursorMz, "Text")
            out.scan(iScan).precursorMz = renameStructField(out.scan(iScan).precursorMz, "Text", "value");
        end
    end
else
    out.scan = formatScan();
end

% Format top level out struct
if isfield(out, "Text")
    out = rmfield(out, "Text");
end
[C, ia] = intersect(["scan" "mzXML" "index"], fieldnames(out));
out = orderfields(out, C(ia));

% Impose NVP constraints/filters
if isfield(nvpArgs, "Levels")
    keepIdx = ismember([out.scan.msLevel], nvpArgs.Levels);
    out.scan(~keepIdx) = [];
end

if isfield(nvpArgs, "TimeRange")
    nvpArgs.TimeRange = sort(nvpArgs.TimeRange);
    scanTimes = arrayfun(@(x) sscanf(x, 'PT%f'), [out.scan.retentionTime]');
    
    if nvpArgs.TimeRange(1) < min(scanTimes) || nvpArgs.TimeRange(2) > max(scanTimes)
        error(message('bioinfo:mzxmlread:InvalidTimeRange',...
            num2str(min(scanTimes)), num2str(max(scanTimes))));
    end
    
    keepIdx = scanTimes > nvpArgs.TimeRange(1) & scanTimes < nvpArgs.TimeRange(2);
    out.scan(~keepIdx) = [];
end

if nvpArgs.Verbose
    disp("DONE")
end
end


%% Helper functions
function outputStruct = addExtraFields(inputStruct, neededFieldNames)
    % add extra fields to the inputStruct, with value []
    outputStruct = inputStruct;
    inputFieldNames = fieldnames(inputStruct);
    extraFieldNames = setdiff(neededFieldNames, inputFieldNames);
    for iF = 1:numel(extraFieldNames)
        outputStruct(1).(extraFieldNames{iF}) = [];
    end
end

function topLevelFields = getTopLevelFieldNamess(scans)
    % return all top-level fields names in all scans (including scans within scans)
    if ~isstruct(scans)
        topLevelFields = {};
        return
    end

    topLevelFields = fieldnames(scans);

    if isfield(scans, "scan") &&  all(cellfun(@(x) isstruct(x) || ismissing(x), {scans.scan}))
        % if scan field value is missing or is a struct, it is not an
        % attribute and can be removed
        topLevelFields = setdiff(topLevelFields, {'scan'});
    end
    
    if isfield(scans, "scan") && any(cellfun(@isstruct, {scans.scan}))
        % get top-level field names of scans within scans
        subScanFields = arrayfun(@(x) getTopLevelFieldNamess(x.scan), scans, "UniformOutput",false);
        subScanFields = [subScanFields{:}];
        topLevelFields = union(subScanFields, topLevelFields);
    end
end


function parsedScan = parseScan(scan, scanIndices)
    % format and flatten scan structs (e.g. pull scans-within-scans to top level)
    parsedScan = [];
    subScan = [];

    if isfield(scan, "scan") && isstruct(scan.scan)
        % remove child element scans, they will be parsed recursively below
        subScan = scan.scan;
        scan = rmfield(scan, "scan");
    end

    if isfield(scan, "num") && ismember(scan.num, scanIndices)
        % format scans within "scanIndices"
        parsedScan = formatScan(scan);
    end

    for iS = 1:numel(subScan)
        % recursively call parseScan for subScans and merge scans
        parsedScan = mergeScans(parsedScan, parseScan(subScan(iS), scanIndices));
    end
end

function formattedScan = formatScan(scan)
    % mzxmlread scans must have the following fields. Additional attributes
    % (fields) will be added to the formattedScan
    %
    % for "precursorMz" and "peaks", the "Text" field in the substruct is
    % the inner text in the xml element. This name is changed to "value"
    % after parseScan
    templateScan = struct('num',[],'msLevel',[],'peaksCount',[],...
        'polarity',[],'scanType',[],'centroided',[],'deisotoped',[],...
        'chargeDeconvoluted',[],'retentionTime',[],'ionisationEnergy',[],...
        'collisionEnergy',[],'collisionGas',[],'collisionGasPressure',[],...
        'startMz',[],'endMz',[],'lowMz',[],'highMz',[],'basePeakMz',[],...
        'basePeakIntensity',[],'totIonCurrent',[],...
        'scanOrigin',struct('parentFileID',[],'num',[]),...
        'precursorMz',struct('precursorScanNum',[],'precursorIntensity',[],'precursorCharge',[],'windowWideness',[],'Text',[]),...
        'maldi',[],'peaks',struct('precision',[],'byteOrder',[],'pairOrder',[],'Text',[]),...
        'nameValue',struct('name',[],'value',[],'type',[]),...
        'comment',[]);
    
    if nargin < 1
        % create empty (0x0) scan struct with template fields
        formattedScan = templateScan;
        formattedScan = formattedScan([]);
        return
    end
    
    formattedScan = addValues2StructTemplate(scan, templateScan);
end

function outputStruct = addValues2StructTemplate(inputStruct, structTemplate)
    % Adds values from inputStruct to structTemplate. Include fields that
    % are not in structTemplate
    arguments
        inputStruct
        structTemplate struct = struct()
    end
    outputStruct    = structTemplate;
    scanFieldNames  = string(fieldnames(inputStruct));
    fieldIsStruct   = structfun(@isstruct, inputStruct);

    for iF = 1:numel(scanFieldNames)
        if fieldIsStruct(iF)
            % if field is struct, recursively call addValues2StructTemplate
            % if structTemplate has the field, use this value as the
            % template. Otherwise, no template (empty struct)
            if isfield(structTemplate, scanFieldNames(iF))
                outputStruct.(scanFieldNames(iF)) = addValues2StructTemplate(inputStruct.(scanFieldNames(iF)), outputStruct.(scanFieldNames(iF)));
            else
                outputStruct.(scanFieldNames(iF)) = addValues2StructTemplate(inputStruct.(scanFieldNames(iF)));
            end
        elseif ~ismissing(inputStruct.(scanFieldNames(iF)))
            % readstruct default value is missing when an attribute is
            % present in one struct, but not present in another. For fields
            % in structTemplate, leave as default structTemplate value
            % ([]). Otherwise the value will be set to [] in mergeScans if
            % the value is top-level
            outputStruct.(scanFieldNames(iF)) = inputStruct.(scanFieldNames(iF));
        end
    end
end

function mzpeaks = processPeaks(peaks, precision)
    % top level shared variables/objects
    
    % The number of elements in peaks must be a multiple of 8 as we are
    % converting it to single or double precision floating point numbers.
    peaks = convertStringsToChars(peaks);
    if mod(numel(peaks), 8) ~= 0
        error(message('bioinfo:mzxmlread:InvalidNumberPeaks', numel(peaks)));                    
    end
    
    switch precision
        case 32
            mzpeaks = typecast(peaks, 'single');
        case 64
            mzpeaks = typecast(peaks, 'double');
        otherwise
            warning(message('bioinfo:mzxmlread:invalidPrecis', precision));
            mzpeaks = typecast(peaks, 'single');
    end
    
    [~, ~, endian] = computer;
    
    if strcmp(endian, 'L')
        mzpeaks = swapbytes(mzpeaks);
    end
end %processPeaks

function outputStruct = renameStructField(inputStruct, oldFieldName, newFieldName)
    outputStruct = inputStruct;
    
    % Get field order
    fieldOrder = replace(fieldnames(inputStruct), oldFieldName, newFieldName);
    
    % Rename field
    [outputStruct.(newFieldName)] = inputStruct.(oldFieldName);
    outputStruct = rmfield(outputStruct, oldFieldName);
    
    % Re-order fields
    outputStruct = orderfields(outputStruct, fieldOrder);
end

function mergedScan = mergeScans(scan1, scan2)
    % mergeScans finds the set difference between s1 and s2, adds fields with
    % empty value to each struct, then concatenates the structs
    
    if isstruct(scan1) && isstruct(scan2)
        scan1Fieldnames = fieldnames(scan1);
        scan2Fieldnames = fieldnames(scan2);

        [~, scan1UniqueIdx, scan2UniqueIdx]  = setxor(scan1Fieldnames, scan2Fieldnames);

        for iS1 = scan1UniqueIdx(:)'
            scan2(1).(scan1Fieldnames{iS1}) = [];
        end

        for iS2 = scan2UniqueIdx(:)'
            scan1(1).(scan2Fieldnames{iS2}) = [];
        end
    end
    mergedScan = [scan1; scan2];
end