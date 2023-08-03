function [mzroi,MSroi,roicell]=ROIpeaks(peaks,thresh,mzerror,minroi,nrows,time)
% function [mzroi,MSroi,roicell]=ROIpeaksnew(peaks,thresh,mzerror,minroi,nrows,time)
% 
% This function allows building a MS data matrix from structure variable peaks 
% selecting only the regions of interest (roi). These regions are defined acording to 
% the following input parameters: thresh, mzeror and minroi 
% 
% INPUT
% % peaks is the cell variable containing MS measurements with as many cells 
% as MS spectra/retention times (nrows). In every cell, mz and MS intensities 
% are given for every spectrum (every cell/MS spectrum has different lengths)
% % thresh is a parameter to filter significative MS intensities, 
% i.e thresh = 0.1% max measured intensity (max(max(MSroi))
% % mzerror is a parameter to define the width of mz experimental values in peaks
% to be considered coming from the same theoreical mz value 
% % minroi minimum number of elution times to be considered in a roi (e.g. 3)
% % nrows number of cells/rows/spectra of the variable peaks to be processed
% % time elution (retention) times corresponding to cells/rows/spectra 
%
% OUTPUT
% % mzroi, finally selected mz values (in the same roicell) 
% % MSroi, new arrenged MS spectra data matrix of dimensions 
%   (nr.of MS spectra), nr.of mzroi values )
% % roicell contents of selected ROI arranged in a cell data set, 
% in the output they are finally ordered according to mzroi values
% roicell{:,1}= individual measured mz values in every ROI (one 
% roicell(:,2}= elution times where couples mz,MSI are measured
% roicell{:,3}= measure MSI intensities at time, mz
% roicell{:,4}= consiedered row of peaks file (elution time);
% roicell{:,5}= mzroi, final mz balue of the considered ROI (mean of all mz values 
% included in ROI 
%
% e.g. mzroi,MSroi,roicell]=ROIpeaksnew(peaks,1000,0.01,10,1899,time);
% where thresh=1000, mzerror=0.01 and minroi=10
% background in MSroi is MSroi=randn(nrows,nmzroi).*0.3*thresh;
%%
% [mzroi,MSroi,roicell]=ROIpeaks(peaks,thresh,mzerror,minroi,nrows,time)

% Adjustments: Addition of random noise removed

mzroi=[];
MSroi=[];
roicell{1,1}=[];
roicell{1,2}=[];
roicell{1,3}=[];
roicell{1,4}=[];
nmzroi=1;

% looking for mzroi values  


for irow=1:nrows
    A=cell2mat(peaks(irow,1));
    A=double(A);
    ipeak=find(A(:,2)>thresh);
    if isfinite(ipeak)
        mz=A(ipeak,1);
        MS=A(ipeak,2);
        if irow==1,mzroi=mz(1);end
        
        
        nmz=size(mz);
        
        for i=1:nmz
           
            %ieq=find(abs(mzroi-mz(i))<=mzerror/2);
            ieq=find(abs(mzroi-mz(i))<=mzerror);
                     
            if isfinite(ieq)
                ieq=ieq(1); 
                roicell{ieq,1}=[roicell{ieq,1},mz(i)];
                roicell{ieq,2}=[roicell{ieq,2},time(irow)];
                roicell{ieq,3}=[roicell{ieq,3},MS(i)];
                roicell{ieq,4}=[roicell{ieq,4},irow];
                roicell{ieq,5}=mean(roicell{ieq,1});
                mzroi(ieq)=roicell{ieq,5};
                
            else
                
                nmzroi=nmzroi+1;
                roicell{nmzroi,1}=mz(i);
                roicell{nmzroi,2}=time(irow);
                roicell{nmzroi,3}=MS(i);
                roicell{nmzroi,4}=irow;
                roicell{nmzroi,5}=mz(i);
                mzroi(nmzroi)=mz(i);
                
            end
            
        end
        
    end
    
    
end

% sort mzroi values
[mzroi,isort]=sort(mzroi);

for i=1:nmzroi,for j=1:5,roicellsort{i,j}=roicell{isort(i),j};end,end
roicell=roicellsort;

% Now, filter those having a minimum number of elution times (minroi)
% and a maximum value higher than thresh

for i=1:nmzroi 
    if isempty(roicell{i,1}),roicell{i,1}=0;end
    numberroi(i)=length(roicell{i,1});
end

for i=1:nmzroi
    if isempty(roicell{i,3}),roicell{i,3}=0;end
    maxroi(i)=max(roicell{i,3});
end

iroi=find(numberroi>minroi & maxroi>thresh);

mzroi=mzroi(iroi);
nmzroi=length(mzroi);
roicell=roicell(iroi,:);

% Evaluation of MS values from roicell{nmzroi,3}
% Now evaluating MS matrix only for thes mzroi values
% Defining first the backgound
MSroi=zeros(nrows,nmzroi);

for i=1:nmzroi
    nval=length(roicell{i,4});
    for j=1:nval
        irow=roicell{i,4}(j);
        MSI=roicell{i,3}(j);
        MSroi(irow,i)=MSroi(irow,i)+MSI;
    end
    % Fill MSroi empty times in mzroi, specially those with zeros inside 
    % a chromatograhic peak
    y=MSroi(:,i);
    iy=find(y>0);
    if numel(iy)>1
        intertime=[time(iy(1):iy(end))];
        ynew = interp1(time(iy), y(iy), intertime);
        MSroi([iy(1):1:iy(end)],i)=ynew;
    else
        MSroi(:,i)=0;
    end
end

    

    

        