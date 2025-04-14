function [ScanDataMS1,ScanDataMS2] = readmzXML_MSandMS2(DataPath)
%% reads mzXML files and outputs mz, Intensity and time data
% .mzXML files can be 32bit and 64bit encoded and can use zlib compression
% Input: DataPath: data path to .mzXML file
% Outputs:
% ScanDataMS1: struct containing peak data (cell array, with two column
% matrix mass and intensity), retention time and scan polarity
% ScanDataMS2: struct containing MS/MS data. Peak data (cell array, with two column
% matrix mass and intensity), retention time, scan polarity, precursor
% mass, fragmentation energy and fragmentation type

arguments
    DataPath            (1,1) string
end

ScanDataMS1 = struct('profileDataMS1',[],...
    'timeDataMS1',[],...
    'polarityMS1',[]);

ScanDataMS2 = struct('profileDataMS2',[],...
    'timeDataMS2',[],...
    'polarityMS2',[],...
    'precursorMass',[],...
    'fragmentationEnergy',[],...
    'fragmentationType',[]);

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
scanData = cell(scanCountValue,1);

% Loop through all scan elements
for i = 1:(scanNodes.getLength)
    scanNode = scanNodes.item(i-1);

    % Extract msLevel, polarity, and retentionTime attributes from the scan element
    msLevels(i) = str2double(scanNode.getAttribute('msLevel'));
    if msLevels(i) > 2
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
    if msLevels(i) > 2
        scanData{i} = 0;
    else
        if compression == true
            intlist = decodeCompressed(PeaksString,precision,endian);
            scanData{i} = reshape(intlist,2,[])';
        else
            intlist = decodeUncompressed(PeaksString,precision,endian);
            scanData{i} = reshape(intlist,2,[])';
        end
    end
end
% remove empty rows
idx = cellfun(@isempty,scanData);
msLevels(idx) = [];
CollisionEnergy(idx) = [];
FragMethod(idx) = [];
scanData(idx) = [];
polarity(idx) = [];
PrecursorMass(idx) = [];
retentionTime(idx) = [];

%convert times to double and round to first decimal
retentionTime = cellfun(@(x) sscanf(x,'PT %f'), retentionTime);

%convert precursormass to double
PrecursorMass = str2double(PrecursorMass);

%store in output structs
idMS1 = msLevels == 1;
ScanDataMS1.profileDataMS1 = scanData(idMS1);
ScanDataMS1.timeDataMS1 = retentionTime(idMS1);
ScanDataMS1.polarityMS1 = polarity(idMS1);

idMS2 = msLevels == 2;
ScanDataMS2.profileDataMS2 = scanData(idMS2);
ScanDataMS2.timeDataMS2 = retentionTime(idMS2);
ScanDataMS2.polarityMS2 = polarity(idMS2);
ScanDataMS2.precursorMass = PrecursorMass(idMS2);
ScanDataMS2.fragmentationEnergy = CollisionEnergy(idMS2);
ScanDataMS2.fragmentationType = FragMethod(idMS2);

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