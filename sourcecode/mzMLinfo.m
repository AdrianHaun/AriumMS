function [FileInfo,retentionTimes,TIC,BPC,polarity] = mzMLinfo(dataPath)
%Reads basic file information from mzML files, only mzML files without
%numPress compression are supported.
%Output Fileinfo is a struct with fields: FileName, FileModDate, FileSize,
%NumberOf Scans, StartTime, EndTime. Start and End times are in seconds

arguments
    dataPath  (1,1) string
end

FileInfo =  struct('numberOfScansMS1',[],...
                'numberOfScansMSn',[],...
                'startTime',[],...
                'endTime',[],...
                'scanFrequenceMS1',[],...
                'scanFrequenceMS2',[]);

doc = xmlread(dataPath);
spectrumList = doc.getElementsByTagName('spectrumList').item(0);
% Check if the spectrum element exists
if ~isempty(spectrumList)
    % Get the value of the scanCount attribute
    scanCountValue = str2double(spectrumList.getAttribute('count'));
    spectrumNodes = doc.getElementsByTagName('spectrum');
    firstNode = spectrumNodes.item(0);
    %get time unit and msPolarity
    scanElement=firstNode.getElementsByTagName('scanList');
    scanElement=scanElement.item(0).getElementsByTagName('scan');
    scanElement=scanElement.item(0).getElementsByTagName('cvParam');
    retentionTimeUnit = string(scanElement.item(0).getAttribute('unitName'));
else
    error("Corrupt or empty file.")
end

msLevels = zeros(scanCountValue,1);
retentionTimes = zeros(scanCountValue,1);
TIC = zeros(scanCountValue,1);
BPC = zeros(scanCountValue,1);
polarity = strings(scanCountValue,1);

% Iterate through each 'spectrum' element
for iScan = 0:spectrumNodes.getLength - 1
    spectrumElement = spectrumNodes.item(iScan);
    % Extract msLevel and retentionTime attributes from the scan element
    ScanInfos = spectrumElement.getElementsByTagName('cvParam');
    for item = 0:ScanInfos.getLength -1
        attribute = string(ScanInfos.item(item).getAttribute('name'));
        switch attribute
            case "ms level"
                msLevels(iScan+1) = str2double(ScanInfos.item(item).getAttribute('value'));
            case {"positive scan","negative scan"}
                polarity(iScan+1) = attribute;
            case "base peak intensity"
                BPC(iScan+1) = str2double(ScanInfos.item(item).getAttribute('value'));
            case "total ion current"
                TIC(iScan+1) = str2double(ScanInfos.item(item).getAttribute('value'));
            otherwise
                continue
        end
            scanElement = spectrumElement.getElementsByTagName('scan');
            scanElement = scanElement.item(0).getElementsByTagName('cvParam');
            retentionTimes(iScan+1) = str2double(scanElement.item(0).getAttribute('value'));
    end
end
% convert retentionTimes to seconds
switch retentionTimeUnit
    case "minute"
        retentionTimes = retentionTimes*60;
end
MS1ID = msLevels == 1;
TIC(~MS1ID) = [];
BPC(~MS1ID) = [];
retentionTimes(~MS1ID) = [];
polarity(~MS1ID) = [];
polarity(polarity == "positive scan") = "+";
polarity(polarity == "negative scan") = "-";


% Remove entries where TIC is zero
MS1ID = TIC ~= 0;
TIC = TIC(MS1ID);
BPC = BPC(MS1ID);
retentionTimes = retentionTimes(MS1ID);
polarity = polarity(MS1ID);

FileInfo.numberOfScansMS1 = sum(MS1ID);
FileInfo.numberOfScansMSn = sum(~MS1ID);
FileInfo.startTime = retentionTimes(1);
FileInfo.endTime = retentionTimes(end);
end
