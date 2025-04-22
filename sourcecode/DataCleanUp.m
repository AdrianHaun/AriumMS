function [Peak,Time] = DataCleanUp(Peak,Time)
%Remove m/z with Intensity < 100 counts
%   Removes all m/z values with intensity below thresh.

%% Clean Data
parfor j = 1:height(Peak)
    data = Peak{j,1};
    [~,edges] = histcounts(data(:,2));
    idx = data(:,2) <= edges(2)/2;
    data(idx,:)=[];
    Peak{j,1}=data;
end
idx=cellfun(@isempty, Peak);
Peak(idx) = [];
Time(idx) = [];
end
