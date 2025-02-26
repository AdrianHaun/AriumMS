function ScanDataMS1 = readmzCDF(dataPath)

arguments
    dataPath            (1,1) string
end

ScanDataMS1 = struct('profileDataMS1',[],...
    'timeDataMS1',[],...
    'polarityMS1',[]);


mzCDFStruct = mzcdfread(dataPath,'Verbose',false);

[peakList, retentionTimes] = mzcdf2peaks(mzCDFStruct);


ScanDataMS1.peakDataMS1 = peakList;
ScanDataMS1.timeDataMS1 = retentionTimes;
