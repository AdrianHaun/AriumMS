function [FileInfo,retentionTimes,TIC,BPC,polarity] = mzXMLinfo(dataPath)
%Reads basic file information from mzML files, only mzML files without
%numPress compression are supported.
%Output Fileinfo is a struct with fields: FileName, FileModDate, FileSize,
%NumberOf Scans, StartTime, EndTime. Start and End times are in seconds

arguments
    dataPath  (1,1) string {mustBeFile}
end

FileInfo =  struct('numberOfScansMS1',[],...
                'numberOfScansMSn',[],...
                'startTime',[],...
                'endTime',[],...
                'scanFrequenceMS1',[],...
                'scanFrequenceMS2',[]);

%check for empty datafile
try doc = xmlread(dataPath);
    
catch exception
    error('mzXMLinfo:emptyFile', 'Data file is empty or corrupt')   
end

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

msLevels = zeros(scanCountValue,1);
retentionTimes = strings(scanCountValue,1);
TIC = zeros(scanCountValue,1);
BPC = zeros(scanCountValue,1);
polarity = strings(scanCountValue,1);

% Loop through all scan elements
for iScan = 1:(scanNodes.getLength)
    scanNode = scanNodes.item(iScan-1);
    % Extract msLevel, polarity, and retentionTime attributes from the scan element
    msLevels(iScan) = str2double(scanNode.getAttribute('msLevel'));
    retentionTimes(iScan) = scanNode.getAttribute('retentionTime');
    polarity(iScan) = scanNode.getAttribute('polarity');
    TIC(iScan) = str2double(scanNode.getAttribute('totIonCurrent'));
    BPC(iScan) = str2double(scanNode.getAttribute('basePeakIntensity'));
end
%convert times to double and round to first decimal
retentionTimes = cellfun(@(x) sscanf(x,'PT %f'), retentionTimes);

MS1ID = msLevels == 1;
TIC(~MS1ID) = [];
BPC(~MS1ID) = [];
retentionTimes(~MS1ID) = [];
polarity(~MS1ID) = [];
FileInfo.numberOfScansMS1 = sum(MS1ID);
FileInfo.numberOfScansMSn = sum(~MS1ID);
FileInfo.startTime = retentionTimes(1);
FileInfo.endTime = retentionTimes(end);
end
