function [FileInfo,RetentionTimes,TIC,BPC] = mzXMLinfo(DataPath)
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
msRunNode = doc.getElementsByTagName('msRun').item(0);
% Check if the msRun element exists
if ~isempty(msRunNode)
    % Get the value of the scanCount attribute
    scanCountValue = str2double(msRunNode.getAttribute('scanCount'));
    % Find the scan element
    scanNodes = doc.getElementsByTagName('scan');
    FirstScan = scanNodes.item(0);
    polarity = string(FirstScan.getAttribute('polarity'));
else
    error("Corrupt or empty file.")
end
switch polarity
    case "+"
        FileInfo.Polarity = "positive";
    case "-"
        FileInfo.Polarity = "negative";
end

msLevels = zeros(scanCountValue,1);
RetentionTimes = strings(scanCountValue,1);
TIC = zeros(scanCountValue,1);
BPC = zeros(scanCountValue,1);

% Loop through all scan elements
for i = 1:(scanNodes.getLength)
    scanNode = scanNodes.item(i-1);
    % Extract msLevel, polarity, and retentionTime attributes from the scan element
    msLevels(i) = str2double(scanNode.getAttribute('msLevel'));
    RetentionTimes(i) = scanNode.getAttribute('retentionTime');
    TIC(i) = str2double(scanNode.getAttribute('totIonCurrent'));
    BPC(i) = str2double(scanNode.getAttribute('basePeakIntensity'));
end
%convert times to double and round to first decimal
RetentionTimes = cellfun(@(x) sscanf(x,'PT %f'), RetentionTimes);

MS1ID = msLevels == 1;
TIC(~MS1ID)=[];
BPC(~MS1ID)=[];
RetentionTimes(~MS1ID)=[];
FileInfo.NumberOfScansMS1 = sum(MS1ID);
FileInfo.NumberOfScansMSn = sum(~MS1ID);
FileInfo.StartTime = RetentionTimes(1);
FileInfo.EndTime = RetentionTimes(end);
end
