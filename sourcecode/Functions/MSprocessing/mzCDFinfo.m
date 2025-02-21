function [FileInfo,RetentionTimes,TIC,BPC,polarityCells] = mzCDFinfo(dataPath)

arguments
    dataPath  (1,1) string
end

FileInfo = struct('NumberOfScansMS1',[],...
    'NumberOfScansMSn',[],...
    'StartTime',[],...
    'EndTime',[]);


fileInfos = mzcdfinfo(dataPath);

FileInfo.NumberOfScansMS1 = fileInfos.NumberOfScans;
FileInfo.NumberOfScansMSn = NaN;
FileInfo.StartTime = fileInfos.StartTime;
FileInfo.EndTime = fileInfos.EndTime;

mzCDFStruct = mzcdfread(dataPath,'Verbose',false);

polarityCells = repmat({mzCDFStruct.test_ionization_polarity},fileInfos.NumberOfScans,1);
[Peaklist, RetentionTimes] = mzcdf2peaks(mzCDFStruct);

TIC = cellfun(@(x) sum(x(:,2)),Peaklist);
BPC = cellfun(@(x) max(x(:,2)),Peaklist);