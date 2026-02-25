function [scanData,retentionTime,varargout] = readmzML(dataPath,options)
%% reads mzML files and outputs mz, Intensity and time data
% .mzML files can be 32bit and 64bit encoded can use zlib compression, numpress compression is not yet
% supported
% Input: DataPath: data path to .mzML file
% optional inputs:
% MSLevel: specifies the MS level to to extract, a msLevel > 1 enables additional outputs. default is 1
% Outputs:
% Peaks: cell array each containing a two column matrix, with mz values and corresponding intensity
% times: vector with scan times corresponding to each cell of Peaks.
% optional outputs: Scan Polarity, Precursor Mass, Fragmentation Energy, Fragmentation Type
arguments
    dataPath  (1,1) string
    options.MSLevel (1,1) {mustBeInteger} = 1
end

doc = xmlread(dataPath);
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
    intBinary = firstMZ.item(1);
    intBinary = intBinary.getElementsByTagName('cvParam');
    intPrecision = string(intBinary.item(0).getAttribute('name'));
    intEncoding = string(intBinary.item(1).getAttribute('name'));
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
for iSpectrum = 0:spectrumNodes.getLength - 1
    spectrumElement = spectrumNodes.item(iSpectrum);
    % Extract msLevel and retentionTime attributes from the scan element
    ScanInfos = spectrumElement.getElementsByTagName('cvParam');
    polarity(iSpectrum+1) = string(ScanInfos.item(2).getAttribute('name'));
    % find ms level value
    maxCount = ScanInfos.getLength-1;
    counter = 0;
    item = "";
    while item ~= "ms level" && counter <= maxCount
        item = string(ScanInfos.item(counter).getAttribute('name'));
        counter = counter + 1;
    end
    msLevels(iSpectrum+1) = str2double(ScanInfos.item(counter-1).getAttribute('value'));
    if msLevels(iSpectrum+1) ~= options.MSLevel
        continue
    else
        scanElement = spectrumElement.getElementsByTagName('scan');
        scanElement = scanElement.item(0).getElementsByTagName('cvParam');
        retentionTime(iSpectrum+1) = str2double(scanElement.item(0).getAttribute('value'));
        % Extract mz and Int binaries
        Binaries = spectrumElement.getElementsByTagName('binaryDataArray');
        mzBinary = Binaries.item(0);
        mzBinaryStrings(iSpectrum+1) = mzBinary.getTextContent;
        intBinary = Binaries.item(1);
        intBinaryStrings(iSpectrum+1) = intBinary.getTextContent;
        if msLevels(iSpectrum+1) > 1
            precursorElement = spectrumElement.getElementsByTagName('precursorList');
            precursorNode = precursorElement.item(0).getElementsByTagName('selectedIonList');
            precursorNode = precursorNode.item(0).getElementsByTagName('cvParam');
            precursorMass(iSpectrum+1) = precursorNode.item(0).getAttribute('value');
            fragmentationNode = precursorElement.item(0).getElementsByTagName('activation');
            fragmentationNode = fragmentationNode.item(0).getElementsByTagName('cvParam');
            fragMethod(iSpectrum+1) = fragmentationNode.item(0).getAttribute('name');
            collisionEnergy(iSpectrum+1) = str2double(fragmentationNode.item(1).getAttribute('value'));
        end
    end
end

%process scan data
%remove white space from binary strings
intBinaryStrings = strtrim(intBinaryStrings);
mzBinaryStrings = strtrim(mzBinaryStrings);
%remove MSn data if not relevant
idx = msLevels ~= options.MSLevel;
retentionTime(idx) = [];
polarity(idx) = [];
collisionEnergy(idx) = [];
fragMethod(idx) = [];
precursorMass(idx) = [];
mzBinaryStrings(idx) = [];
intBinaryStrings(idx) = [];
% convert retentionTimes to seconds
switch retentionTimeUnit
    case "minute"
        retentionTime = retentionTime*60;
end
scanData = cell(size(mzBinaryStrings,1),1);

for iScan=1:size(mzBinaryStrings,1) %decode binary strings
    %mzBinary
    if mzEncoding == true
        mzbinary = decodeCompressed(mzBinaryStrings(iScan),mzPrecision);
    else
        mzbinary = decodeUncompressed(mzBinaryStrings(iScan),mzPrecision);
    end
    if intEncoding == true
        Intbinary = decodeCompressed(intBinaryStrings(iScan),intPrecision);
    else
        Intbinary = decodeUncompressed(intBinaryStrings(iScan),intPrecision);
    end
    scanData{iScan} = [mzbinary' Intbinary'];
end

% convert polarity
polarity(polarity == "positive scan") = "+";
polarity(polarity == "negative scan") = "-";

% remove empty rows
idx = cellfun(@isempty,scanData);
collisionEnergy(idx) = [];
fragMethod(idx) = [];
scanData(idx) = [];
polarity(idx) = [];
precursorMass(idx) = [];
retentionTime(idx) = [];

varargout{1} = polarity;

if options.MSLevel > 1
    varargout{2} = str2double(precursorMass);
    varargout{3} = collisionEnergy;
    varargout{4} = fragMethod;

else
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = [];
end


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
