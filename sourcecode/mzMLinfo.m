function [FileInfo,RetentionTimes,TIC,BPC] = mzMLinfo(DataPath)
%Reads basic file information from mzML files, only mzML files without
%numPress compression are supported.
%Output Fileinfo is a struct with fields: FileName, FileModDate, FileSize,
%NumberOf Scans, StartTime, EndTime. Start and End times are in seconds

arguments
    DataPath  (1,1) string
end

FileInfo = struct('Polarity','N/A',...
    'NumberOfScansMS1',[],...
    'NumberOfScansMSn',[],...
    'StartTime',[],...
    'EndTime',[]);

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
else
    error("Corrupt or empty file.")
end

switch polarity
    case "positive scan"
        FileInfo.Polarity = "positive";
    case "negative scan"
        FileInfo.Polarity = "negative";
end

msLevels = zeros(scanCountValue,1);
RetentionTimes = zeros(scanCountValue,1);
TIC = zeros(scanCountValue,1);
BPC = zeros(scanCountValue,1);

% Iterate through each 'spectrum' element
for i = 0:spectrumNodes.getLength - 1
    spectrumElement = spectrumNodes.item(i);
    % Extract msLevel and retentionTime attributes from the scan element
    ScanInfos = spectrumElement.getElementsByTagName('cvParam');
    msLevels(i+1) = str2double(ScanInfos.item(1).getAttribute('value'));
    BPC(i+1) = str2double(ScanInfos.item(5).getAttribute('value'));
    TIC(i+1) = str2double(ScanInfos.item(6).getAttribute('value'));
    scanElement = spectrumElement.getElementsByTagName('scan');
    scanElement = scanElement.item(0).getElementsByTagName('cvParam');
    RetentionTimes(i+1) = str2double(scanElement.item(0).getAttribute('value'));
end
% convert retentionTimes to seconds
switch retentionTimeUnit
    case "minute"
        RetentionTimes = RetentionTimes*60;
end
MS1ID = msLevels == 1;
TIC(~MS1ID)=[];
BPC(~MS1ID)=[];
RetentionTimes(~MS1ID)=[];
FileInfo.NumberOfScansMS1 = sum(MS1ID);
FileInfo.NumberOfScansMSn = sum(~MS1ID);
FileInfo.StartTime = RetentionTimes(1);
FileInfo.EndTime = RetentionTimes(end);
end
