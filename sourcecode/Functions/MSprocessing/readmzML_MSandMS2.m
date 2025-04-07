function [scanDataMS1,scanDataMS2] = readmzML_MSandMS2(dataPath)
%% reads mzML files and outputs mz, Intensity and time data
% .mzML files can be 32bit and 64bit encoded can use zlib compression, numpress compression is not yet
% supported
% Input: DataPath: data path to .mzML file
%
% Outputs:
% ScanDataMS1: struct containing peak data (cell array, with two column
% matrix mass and intensity), retention time and scan polarity
% ScanDataMS2: struct containing MS/MS data. Peak data (cell array, with two column
% matrix mass and intensity), retention time, scan polarity, precursor
% mass, fragmentation energy and fragmentation type

arguments
    dataPath            (1,1) string {mustBeFile}
end

scanDataMS1 = struct('profileDataMS1',[],...
    'timeDataMS1',[],...
    'polarityMS1',[]);

scanDataMS2 = struct('profileDataMS2',[],...
    'timeDataMS2',[],...
    'polarityMS2',[],...
    'precursorMass',[],...
    'fragmentationEnergy',[],...
    'fragmentationType',[]);

doc = xmlread(dataPath);
spectrumList = doc.getElementsByTagName('spectrumList').item(0);
% Check if the spectrum element exists
if isempty(spectrumList)
    error("Corrupt or empty file.")
else
    % Get the value of the scanCount attribute
    scanCountValue = str2double(spectrumList.getAttribute('count'));
    spectrumNodes = doc.getElementsByTagName('spectrum');
    firstNode = spectrumNodes.item(0);
    %get time unit
    scanElement = firstNode.getElementsByTagName('scanList');
    scanElement = scanElement.item(0).getElementsByTagName('scan');
    scanElement = scanElement.item(0).getElementsByTagName('cvParam');
    retentionTimeUnit = string(scanElement.item(0).getAttribute('unitName'));
    %get encoding precision and compression
    firstMZ = firstNode.getElementsByTagName('binaryDataArray');
    mzBinary = firstMZ.item(0);
    mzBinary = mzBinary.getElementsByTagName('cvParam');
    mzPrecision = string(mzBinary.item(0).getAttribute('name'));
    mzEncoding = string(mzBinary.item(1).getAttribute('name'));
    intBinary = firstMZ.item(1);
    intBinary = intBinary.getElementsByTagName('cvParam');
    intPrecision = string(intBinary.item(0).getAttribute('name'));
    intEncoding = string(intBinary.item(1).getAttribute('name'));
end

switch mzPrecision
    case "32-bit float"
        mzPrecision = 'single';
    case "32-bit integer"
        mzPrecision = 'single';
    case "64-bit float"
        mzPrecision = 'double';
    case "64-bit integer"
        mzPrecision = 'double';
end
switch intPrecision
    case "32-bit float"
        intPrecision = 'single';
    case "32-bit integer"
        intPrecision = 'single';
    case "64-bit float"
        intPrecision = 'double';
    case "64-bit integer"
        intPrecision = 'double';
end
switch mzEncoding
    case "no compression"
        mzEncoding = false;
    case "zlib compression"
        mzEncoding = true;
    otherwise
        error("Files use unsupported compression type")
end
switch intEncoding
    case "no compression"
        intEncoding = false;
    case "zlib compression"
        intEncoding = true;
    otherwise
        error("Files use unsupported compression type")
end

% Initialize arrays to store the extracted data
msLevels = zeros(scanCountValue,1);
retentionTime = zeros(scanCountValue,1);
collisionEnergy = zeros(scanCountValue,1);
polarity = strings(scanCountValue,1);
fragMethod = strings(scanCountValue,1);
precursorMass = strings(scanCountValue,1);
mzBinaryStrings = strings(scanCountValue,1);
intBinaryStrings = strings(scanCountValue,1);

% Iterate through each 'spectrum' element
for i = 0:spectrumNodes.getLength - 1
    spectrumElement = spectrumNodes.item(i);
    % Extract msLevel and retentionTime attributes from the scan element
    ScanInfos = spectrumElement.getElementsByTagName('cvParam');
    polarity(i+1) = string(ScanInfos.item(2).getAttribute('name'));
    % find ms level value
    maxCount = ScanInfos.getLength-1;
    counter = 0;
    item = "";
    while item ~= "ms level" && counter <= maxCount
        item = string(ScanInfos.item(counter).getAttribute('name'));
        counter = counter + 1;
    end
    msLevels(i+1) = str2double(ScanInfos.item(counter-1).getAttribute('value'));
    scanElement = spectrumElement.getElementsByTagName('scan');
    scanElement = scanElement.item(0).getElementsByTagName('cvParam');
    retentionTime(i+1) = str2double(scanElement.item(0).getAttribute('value'));
    % Extract mz and Int binaries
    Binaries = spectrumElement.getElementsByTagName('binaryDataArray');
    mzBinary = Binaries.item(0);
    mzBinaryStrings(i+1) = mzBinary.getTextContent;
    intBinary = Binaries.item(1);
    intBinaryStrings(i+1) = intBinary.getTextContent;
    if msLevels(i+1) > 1
        precursorElement = spectrumElement.getElementsByTagName('precursorList');
        precursorNode = precursorElement.item(0).getElementsByTagName('selectedIonList');
        precursorNode = precursorNode.item(0).getElementsByTagName('cvParam');
        precursorMass(i+1) = precursorNode.item(0).getAttribute('value');
        fragmentationNode = precursorElement.item(0).getElementsByTagName('activation');
        fragmentationNode = fragmentationNode.item(0).getElementsByTagName('cvParam');
        fragMethod(i+1) = fragmentationNode.item(0).getAttribute('name');
        collisionEnergy(i+1) = str2double(fragmentationNode.item(1).getAttribute('value'));
    end
end

%process scan data
%remove white space from binary strings
intBinaryStrings = strtrim(intBinaryStrings);
mzBinaryStrings = strtrim(mzBinaryStrings);

scanData = cell(size(mzBinaryStrings,1),1);
for n=1:size(mzBinaryStrings,1) %decode binary strings
    if mzEncoding == true
        mzbinary = decodeCompressed(mzBinaryStrings(n),mzPrecision);
    else
        mzbinary = decodeUncompressed(mzBinaryStrings(n),mzPrecision);
    end
    if intEncoding == true
        Intbinary = decodeCompressed(intBinaryStrings(n),intPrecision);
    else
        Intbinary = decodeUncompressed(intBinaryStrings(n),intPrecision);
    end
    scanData{n} = [mzbinary' Intbinary'];
end

% convert retentionTimes to seconds
switch retentionTimeUnit
    case "minute"
        retentionTime = retentionTime*60;
end
% convert polarity
polarity(polarity == "positive scan") = "+";
polarity(polarity == "negative scan") = "-";

% remove empty rows
idx = cellfun(@isempty,scanData);
msLevels(idx) = [];
collisionEnergy(idx) = [];
fragMethod(idx) = [];
scanData(idx) = [];
polarity(idx) = [];
precursorMass(idx) = [];
retentionTime(idx) = [];

%store in output struct
idMS1 = msLevels == 1;
scanDataMS1.profileDataMS1 = scanData(idMS1);
scanDataMS1.timeDataMS1 = retentionTime(idMS1);
scanDataMS1.polarityMS1 = polarity(idMS1);

idMS2 = msLevels == 2;
scanDataMS2.profileDataMS2 = scanData(idMS2);
scanDataMS2.timeDataMS2 = retentionTime(idMS2);
scanDataMS2.polarityMS2 = polarity(idMS2);
scanDataMS2.precursorMass = str2double(precursorMass(idMS2));
scanDataMS2.fragmentationEnergy = collisionEnergy(idMS2);
scanDataMS2.fragmentationType = fragMethod(idMS2);

%clean data
[scanDataMS1,scanDataMS2] = cleanRawProfileScans(scanDataMS1,scanDataMS2);

end

%helper functions
function out = decodeUncompressed(DataString,precision)
out = typecast(matlab.net.base64decode(DataString),'uint8');
out = typecast(out,precision);
end

function out = decodeCompressed(dataString,precision)
out = typecast(matlab.net.base64decode(dataString),'uint8');
out = zmat(out,0,'zlib');
out = typecast(out,precision);
end
