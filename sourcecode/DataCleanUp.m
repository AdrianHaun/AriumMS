function [Peaklist,Timelist] =DataCleanUp(Peaklist,Timelist)
%Removes noise m/z by binning, and removing the first bin

%% Clean Data
parfor k = 1 : size(Peaklist,1)
    Peak = Peaklist{k,1};
    Time = Timelist{k,1};
    [nrows,~] = size(Peak);
    if nrows == 1
        Peaklist{k,1}=Peak;
    else
        for j = 1:nrows
            data = Peak{j,1};
            [~,edges] = histcounts(data(:,2));
            idx = data(:,2)<= edges(2);
            data(idx,:)=[];
            Peak{j,1}=data;
        end
        idx=cellfun(@isempty, Peak);
        Peak(idx) = [];
        Time(idx) = [];
        Timelist{k,1} = Time;
        Peaklist{k,1} = Peak;
    end
end
end
