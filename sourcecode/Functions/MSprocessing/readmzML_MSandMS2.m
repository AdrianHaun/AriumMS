function [ScanDataMS1,ScanDataMS2] = readmzML_MSandMS2(DataPath)
%% reads mzML files and outputs mz, Intensity and time data
% .mzML files can be 32bit and 64bit encoded can use zlib compression, numpress compression is not yet
% supported
% Input: DataPath: data path to .mzML file

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

doc = xmlread(DataPath);
spectrumList = doc.getElementsByTagName('spectrumList').item(0);
% Check if the spectrum element exists
if ~isempty(spectrumList)
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
    IntBinary = firstMZ.item(1);
    IntBinary = IntBinary.getElementsByTagName('cvParam');
    IntPrecision = string(IntBinary.item(0).getAttribute('name'));
    IntEncoding = string(IntBinary.item(1).getAttribute('name'));
else
    error("Corrupt or empty file.")
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
switch IntPrecision
    case "32-bit float"
        IntPrecision = 'single';
    case "32-bit integer"
        IntPrecision = 'single';
    case "64-bit float"
        IntPrecision = 'double';
    case "64-bit integer"
        IntPrecision = 'double';
end
switch mzEncoding
    case "no compression"
        mzEncoding = false;
    case "zlib compression"
        mzEncoding = true;
    case "MS-Numpress linear prediction compression"
        error("Files that use numpress compression are not supported at this time.")
    case "MS-Numpress linear prediction compression followed by zlib compression"
        error("Files that use numpress compression are not supported at this time.")
end
switch IntEncoding
    case "no compression"
        IntEncoding = false;
    case "zlib compression"
        IntEncoding = true;
    case "MS-Numpress positive integer compression"
        error("Files that use numpress compression are not supported at this time.")
    case "MS-Numpress positive integer compression followed by zlib compression"
        error("Files that use numpress compression are not supported at this time.")
end
% Initialize arrays to store the extracted data
msLevels = zeros(scanCountValue,1);
retentionTime = zeros(scanCountValue,1);
CollisionEnergy = zeros(scanCountValue,1);
polarity = strings(scanCountValue,1);
FragMethod = strings(scanCountValue,1);
PrecursorMass = strings(scanCountValue,1);
mzBinaryStrings = strings(scanCountValue,1);
IntBinaryStrings = strings(scanCountValue,1);

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
    IntBinary = Binaries.item(1);
    IntBinaryStrings(i+1) = IntBinary.getTextContent;
    if msLevels(i+1) > 1
        PrecursorElement = spectrumElement.getElementsByTagName('precursorList');
        PrecursorNode = PrecursorElement.item(0).getElementsByTagName('selectedIonList');
        PrecursorNode = PrecursorNode.item(0).getElementsByTagName('cvParam');
        PrecursorMass(i+1) = PrecursorNode.item(0).getAttribute('value');
        FragmentationNode = PrecursorElement.item(0).getElementsByTagName('activation');
        FragmentationNode = FragmentationNode.item(0).getElementsByTagName('cvParam');
        FragMethod(i+1) = FragmentationNode.item(0).getAttribute('name');
        CollisionEnergy(i+1) = str2double(FragmentationNode.item(1).getAttribute('value'));
    end
end

%process scan data
%remove white space from binary strings
IntBinaryStrings = strtrim(IntBinaryStrings);
mzBinaryStrings = strtrim(mzBinaryStrings);

scanData = cell(size(mzBinaryStrings,1),1);
for n=1:size(mzBinaryStrings,1) %decode binary strings
    if mzEncoding == true
        mzbinary = decodeCompressed(mzBinaryStrings(n),mzPrecision);
    else
        mzbinary = decodeUncompressed(mzBinaryStrings(n),mzPrecision);
    end
    if IntEncoding == true
        Intbinary = decodeCompressed(IntBinaryStrings(n),IntPrecision);
    else
        Intbinary = decodeUncompressed(IntBinaryStrings(n),IntPrecision);
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

%store in output structs
idMS1 = msLevels == 1;
ScanDataMS1.profileDataMS1 = scanData(idMS1);
ScanDataMS1.timeDataMS1 = retentionTime(idMS1);
ScanDataMS1.polarityMS1 = polarity(idMS1);

idMS2 = msLevels == 2;
ScanDataMS2.profileDataMS2 = scanData(idMS2);
ScanDataMS2.timeDataMS2 = retentionTime(idMS2);
ScanDataMS2.polarityMS2 = polarity(idMS2);
ScanDataMS2.precursorMass = str2double(PrecursorMass(idMS2));
ScanDataMS2.fragmentationEnergy = CollisionEnergy(idMS2);
ScanDataMS2.fragmentationType = FragMethod(idMS2);

end
%helper functions

function out = decodeUncompressed(DataString,precision)
out = typecast(matlab.net.base64decode(DataString),'uint8');
out = typecast(out,precision);
end

function out = decodeCompressed(DataString,precision)
out = typecast(matlab.net.base64decode(DataString),'uint8');
out = zmat(out,0,'zlib');
out = typecast(out,precision);
end
