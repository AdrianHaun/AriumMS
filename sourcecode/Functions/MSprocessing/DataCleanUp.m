function msStruct = DataCleanUp(msStruct)
%Determines the most common Intensity values and removes them

%% determine over all intensity bins
allScans = msStruct.profileDataMS1;
ints = vertcat(allScans{:});
[~,edges] = histcounts(ints(:,2));
cutoff = edges(2); 

%% Clean Data
parfor j = 1:height(allScans)
    data = allScans{j,1};
    idx = data(:,2) <= cutoff;
    data(idx,:) = [];
    allScans{j,1} = data;
end
msStruct.profileDataMS1 = allScans;

end
