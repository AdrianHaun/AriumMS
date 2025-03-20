
ProfileData = cell(12,1);
timeList = cell(12,1);
tstart = tic;
parfor n = 1:12
    [ProfileData{n,1},~] = readmzXML(files(n),MSLevel=1)
    [ProfileData{n,1},timeList{n,1}] = DataCleanUp(ProfileData{n,1},timeList{n,1})
end
tend = toc(tstart);

whos ProfileData
disp("Calculation time: "+tend+ " s")

CentroidData = cell(12,1);
tstart = tic;
parfor n = 1:12
    [CentroidData{n,1},timeList{n,1}] = readmzXML(files(n),MSLevel=1)
    CentroidData{n,1} = CentroidScans(CentroidData{n,1});
end
tend = toc(tstart);

whos CentroidData
disp("Calculation time: "+tend+ " s")