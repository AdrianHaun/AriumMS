function [MSroi_end,mzroi_end,time_end]=AutoROI(Peaklist,Timelist,mzerror,mzErrorUnit,minroi,thresh)
%%AutoROI Performs fully automated ROI search and augment.
%   Designd to be used in conjunction with ROIprocess

%preallocate cell arrays
mzlist = cell(length(Peaklist),1);
MSroilist = cell(length(Peaklist),1);
%ROI search for every Sample
parfor d = 1 : length(Peaklist)
    P= Peaklist{d,1};
    T= Timelist{d,1};
    nrows=length(P);
    [mzlist{d,1},MSroilist{d,1},~]=ROIpeaks2(P,thresh,mzerror,mzErrorUnit,minroi,nrows,T);
end
if length(mzlist) == 1 %Skip Augmentation if only one Sample
    MSroi_end=MSroilist{1,1};
    mzroi_end=mzlist{1,1};
    time_end=Timelist{1,1};
else
    for i = 2:size(Peaklist,1)
        [MSroilist{1,1},mzlist{1,1},Timelist{1,1}] = MSroiaug2(MSroilist{1,1},MSroilist{i,1},mzlist{1,1},mzlist{i,1},mzerror,mzErrorUnit,thresh,Timelist{1,1},Timelist{i,1});
    end
    MSroi_end=MSroilist{1,1};
    mzroi_end=mzlist{1,1};
    time_end=Timelist{1,1};
end
end

