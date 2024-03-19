function [ScanData,retentionTime,varargout] = readmzXML(DataPath,options)
%% reads mzXML files and outputs mz, Intensity and time data
% .mzXML files can be 32bit and 64bit encoded can use zlib compression
% Input: Datapath: datapath to .mzML file
% optional inputs:
% MSLevel: specifies the MS level to to extract, a MSLevel > 1 enables additional outputs. default is 1
% Outputs:
% Peaks: cell array each containing a two column matrix, with mz values and corresponding intensity
% times: vector with scan times corresponding to each cell of Peaks.
% optional outputs: Precursor Mass, Fragmentation Energy, Fragmentation Type

arguments
    DataPath (1,1) string
    options.MSLevel (1,1) {mustBeInteger} = 1
end
[~,~,endian] = computer;
doc = xmlread(DataPath);
%get number of scans
%Find the msRun element
msRunNode = doc.getElementsByTagName('msRun').item(0);

% Check if the msRun element exists
if ~isempty(msRunNode)
    % Get the value of the scanCount attribute
    scanCountValue = str2double(msRunNode.getAttribute('scanCount'));
    % Find the scan element
    scanNodes = doc.getElementsByTagName('scan');
else
    error("Corrupt or empty file.")
end

% Initialize arrays to store extracted values
msLevels = zeros(scanCountValue,1);
polarity = strings(scanCountValue,1);
retentionTime = strings(scanCountValue,1);
CollisionEnergy = zeros(scanCountValue,1);
FragMethod = strings(scanCountValue,1);
PrecursorMass = strings(scanCountValue,1);
ScanData = cell(scanCountValue,1);

% Loop through all scan elements
for i = 1:(scanNodes.getLength)
    scanNode = scanNodes.item(i-1);

    % Extract msLevel, polarity, and retentionTime attributes from the scan element
    msLevels(i) = str2double(scanNode.getAttribute('msLevel'));
    if msLevels(i) ~= options.MSLevel
        continue
    end
    polarity(i) = scanNode.getAttribute('polarity');
    retentionTime(i) = scanNode.getAttribute('retentionTime');
    % Extract values from the precoursor element when MS2 scan
    if msLevels(i) > 1
        PrecursorNode = scanNode.getElementsByTagName('precursorMz').item(0);
        CollisionEnergy(i) = str2double(scanNode.getAttribute('collisionEnergy'));
        FragMethod(i) = string(PrecursorNode.getAttribute('activationMethod'));
        PrecursorMass(i) = str2double(PrecursorNode.getTextContent);
    end
    % Extract values from the peaks element
    peaksNode = scanNode.getElementsByTagName('peaks').item(0); % Assuming there's only one peaks element
    precision = str2double(peaksNode.getAttribute('precision'));
    switch precision
        case 32
            precision = 'single';
        case 64
            precision = 'double';
    end
    compressionType = string(peaksNode.getAttribute('compressionType'));
    switch compressionType
        case "none"
            compression = false;
        case "zlib"
            compression = true;
    end
    %get intensity data
    PeaksString = char(peaksNode.getTextContent);
    %decode
    if msLevels(i) ~= options.MSLevel
        ScanData{i}=0;
    else
        if compression == true
            intlist = decodeCompressed(PeaksString,precision,endian);
            ScanData{i} = reshape(intlist,2,[])';
        else
            intlist = decodeUncompressed(PeaksString,precision,endian);
            ScanData{i} = reshape(intlist,2,[])';
        end
    end
end
% remove empty rows
idx = msLevels ~= options.MSLevel;
CollisionEnergy(idx)=[];
FragMethod(idx)=[];
ScanData(idx)=[];
polarity(idx)=[];
PrecursorMass(idx)=[];
retentionTime(idx)=[];
%convert times to double and round to first decimal
retentionTime = cellfun(@(x) sscanf(x,'PT %f'), retentionTime);

if options.MSLevel > 1
    varargout{1} = str2double(PrecursorMass);
    varargout{2} = CollisionEnergy;
    varargout{3} = FragMethod;
else
    varargout{1} = 0;
    varargout{2} = 0;
    varargout{3} = "empty";
end
end

%helper functions
function out = decodeUncompressed(DataString,precision,endian)
out = typecast(matlab.net.base64decode(DataString),'int8');
out = typecast(out,precision);
if strcmp(endian,'L')
    out = swapbytes(out);
end
end

function out = decodeCompressed(DataString,precision,endian)
out = typecast(matlab.net.base64decode(DataString),'int8');
out = zmat(out,0,'zlib');
out = typecast(out,precision);
if strcmp(endian,'L')
    out = swapbytes(out);
end
end