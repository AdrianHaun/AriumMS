function Peaklist =DataCleanUp(Peaklist)
%Remove m/z with Intensity < 100 counts
%   Removes all m/z values with intensity below thresh.

%% Clean Data
parfor k = 1 : size(Peaklist,1)
    Peak=Peaklist{k,1};
    [nrows,~] = size(Peak);
    if nrows == 1
        Peaklist{k,1}=Peak;
    else
        for j = 1:nrows
            data=Peak{j,1};
            idx=data(:,2)< 100;
            data(idx,:)=[];
            Peak{j,1}=data;
        end
        Peaklist{k,1}=Peak;
    end
end
end