function [MSroi_end,mzroi_end]=filterContaminants(MSroi_end,mzroi_end,mzerror,Contaminants)
%% Removes common contaminants from mass list
%find Contaminants
idx=Contaminants<mzroi_end+mzerror/2 & Contaminants>mzroi_end-mzerror/2;
idx=any(idx,1);
%remove contaminant masses from mzroi_end and MSroi_end
mzroi_end(idx)=[];
MSroi_end(:,idx)=[];
end
