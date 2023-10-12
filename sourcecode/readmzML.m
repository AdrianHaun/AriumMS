function [ScanData,retentionTime,varargout] = readmzML(DataPath,options)
%% reads mzML files and outputs mz, Intensity and time data
% .mzML files can be 32bit and 64bit encoded can use zlib compression, numpress compression is not yet
% supported
% Input: Datapath: datapath to .mzML file
% optional inputs:
% MSLevel: specifies the MS level to to extract, a MSLevel > 1 enables additional outputs. default is 1
% Outputs:
% Peaks: cell array each containing a two column matrix, with mz values and corresponding intensity
% times: vector with scan times corresponding to each cell of Peaks.
% optional outputs: Precursor Mass, Fragmentation Energy, Fragmentation
% Type
arguments
    DataPath            (1,1) string
    options.MSLevel     (1,1) {mustBeInteger} = 1
end

doc = xmlread(DataPath);
spectrumList = doc.getElementsByTagName('spectrumList').item(0);
% Check if the spectrum element exists
if ~isempty(spectrumList)
    % Get the value of the scanCount attribute
    scanCountValue = str2double(spectrumList.getAttribute('count'));
    spectrumNodes = doc.getElementsByTagName('spectrum');
    firstNode = spectrumNodes.item(0);
    %get time unit and msPolarity
    ScanInfos = firstNode.getElementsByTagName('cvParam');
    polarity = string(ScanInfos.item(2).getAttribute('name'));
    scanElement=firstNode.getElementsByTagName('scanList');
    scanElement=scanElement.item(0).getElementsByTagName('scan');
    scanElement=scanElement.item(0).getElementsByTagName('cvParam');
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
%
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
FragMethod = strings(scanCountValue,1);
PrecursorMass = strings(scanCountValue,1);
mzBinaryStrings = strings(scanCountValue,1);
IntBinaryStrings = strings(scanCountValue,1);

% Iterate through each 'spectrum' element
for i = 0:spectrumNodes.getLength - 1
    spectrumElement = spectrumNodes.item(i);
    % Extract msLevel and retentionTime attributes from the scan element
    ScanInfos = spectrumElement.getElementsByTagName('cvParam');
    msLevels(i+1) = str2double(ScanInfos.item(1).getAttribute('value'));
    if msLevels(i+1) ~= options.MSLevel
        continue
    else
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
end
%remove whitespace from binary strings
IntBinaryStrings = strtrim(IntBinaryStrings);
mzBinaryStrings = strtrim(mzBinaryStrings);
%remove MSn data if not relevant
idx = msLevels ~= options.MSLevel;
retentionTime(idx)=[];
CollisionEnergy(idx)=[];
FragMethod(idx)=[];
PrecursorMass(idx)=[];
mzBinaryStrings(idx)=[];
IntBinaryStrings(idx)=[];
% convert retentionTimes to seconds
switch retentionTimeUnit
    case "minute"
        retentionTime = retentionTime*60;
end
ScanData = cell(size(mzBinaryStrings,1),1);
for n=1:size(mzBinaryStrings,1)
    %mzBinary
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
    ScanData{n} = [mzbinary' Intbinary'];
end

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

function out = decodeUncompressed(DataString,precision)
out = typecast(matlab.net.base64decode(DataString),'uint8');
out = typecast(out,precision);
end

function out = decodeCompressed(DataString,precision)
out = typecast(matlab.net.base64decode(DataString),'uint8');
out = zmat(out,0,'zlib');
out = typecast(out,precision);
end
