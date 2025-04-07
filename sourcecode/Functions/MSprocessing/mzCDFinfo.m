function [FileInfo,RetentionTimes,TIC,BPC,polarityCells] = mzCDFinfo(dataPath)

arguments
    dataPath  (1,1) string
end

FileInfo =  struct('numberOfScansMS1',[],...
                'numberOfScansMSn',[],...
                'startTime',[],...
                'endTime',[],...
                'scanFrequenceMS1',[],...
                'scanFrequenceMS2',[]);


fileInfos = mzcdfinfo(dataPath);

FileInfo.numberOfScansMS1 = fileInfos.NumberOfScans;
FileInfo.numberOfScansMSn = NaN;
FileInfo.startTime = fileInfos.StartTime;
FileInfo.endTime = fileInfos.EndTime;

mzCDFStruct = mzcdfread(dataPath,'Verbose',false);

polarityCells = repmat({mzCDFStruct.test_ionization_polarity},fileInfos.NumberOfScans,1);
[Peaklist, RetentionTimes] = mzcdf2peaks(mzCDFStruct);

TIC = cellfun(@(x) sum(x(:,2)),Peaklist);
BPC = cellfun(@(x) max(x(:,2)),Peaklist);