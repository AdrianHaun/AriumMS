function [ScanDataMS1,ScanDataMS2] = readmzCDF(dataPath)

arguments
    dataPath    (1,1) string
end

ScanDataMS1 = struct('profileDataMS1',[],...
    'timeDataMS1',[],...
    'polarityMS1',[]);


mzCDFStruct = mzcdfread(dataPath,'Verbose',false);

[peakList, retentionTimes] = mzcdf2peaks(mzCDFStruct);

%remove possible empty scans
idx = cellfun(@isempty,peakList);
peakList(idx) = [];
retentionTimes(idx) = [];

ScanDataMS1.profileDataMS1 = peakList;
ScanDataMS1.timeDataMS1 = retentionTimes;

%output dummy MS2 data struct for compatibility with other MS file
%processing functions
ScanDataMS2 = struct('profileDataMS2',[],...
    'timeDataMS2',[],...
    'polarityMS2',[],...
    'precursorMass',[],...
    'fragmentationEnergy',[],...
    'fragmentationType',[]);
