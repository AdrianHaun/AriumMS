function msStruct = DataCleanUpMS2(msStruct,cutoff)
%Removes masses with intensity less than cutoff value

allScans = msStruct.centroidDataMS2;

%% Clean Data
parfor j = 1:height(allScans)
    data = allScans{j,1};
    %rescale intensity
    data(:,2) = data(:,2)./max(data(:,2));
    idx = data(:,2) < cutoff;
    data(idx,:) = [];
    allScans{j,1} = data;
end
msStruct.centroidDataMS2 = allScans;

end
