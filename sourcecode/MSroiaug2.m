function [MSroi_aug,mzroi_aug,time_aug] = MSroiaug2(MSroi1,MSroi2,mzroi1,mzroi2,mzerror,mzErrorUnit,thresh,time1,time2)

% INPUT
% MSroi1 and MSroi2, the two MS matrices to augment
% mzroi1,mzroi2, the two mzroi values consoidered in each matrix
% mzerror is a parameter to define the diff between two mazroi values
% to be considerd the same (+- mzerror)
% tresh, threshold value to be considered in the augmented MSroi_aug and mzroi_aug
% time1 and time2 are the two elution time measurements
%
% OUTPUT
% MSroi_aug is the column-wise augmented MS matrix,
% mzroi_aug are the new mzroi values finally considered

% Adjustment: replaced the addition of random noise with 0
% added relative mz error [ppm] support

[nr1,~]=size(MSroi1);
[nr2,~]=size(MSroi2);
mzroi_non=[];
time_aug=[time1;time2];
imz=0;
imzn=0;
% Check for matching mzroi values in first matrix

for i=1:length(mzroi1)
    switch mzErrorUnit
        case "Da"
            j=find(abs(mzroi1(1,i)-mzroi2(1,:))<= mzerror);
        case "ppm"
            j=find(abs(((mzroi1(1,i)-mzroi2(1,:))/mzroi1(1,i))*10^6)<= mzerror);
    end

    if isfinite(j)
        jsize=length(j);
        imz=imz+1;
        if jsize>1
            mzroi_aug(imz)=mean([mzroi1(1,i),mean([mzroi2(1,j(1):j(jsize))])]);
            MS=sum(MSroi2(:,j(1):j(jsize))')';
            MSroi_aug(:,imz)=[MSroi1(:,i);MS];
        else
            mzroi_aug(imz)=mean([mzroi1(1,i),mzroi2(1,j)]);
            MSroi_aug(:,imz)=[MSroi1(:,i);MSroi2(:,j)];
        end
        
    else
        
        % Check now for non-matching mzroi1 values
        
        if max(MSroi1(:,i)) >= thresh
            imz=imz+1;
            mzroi_aug(imz)=mzroi1(1,i);
            MSroi_aug(:,imz)=[MSroi1(:,i);zeros(nr2,1)];
            
        else
            imzn=imzn+1;
            mzroi_non(imzn)=mzroi1(1,i);
            MSroi_non(:,imzn)=[MSroi1(:,i);zeros(nr2,1)];
           
        end
    end
end

for i=1:length(mzroi2)
    j=find(abs(mzroi2(1,i)-mzroi1(1,:))<=mzerror/2);
    
    if isempty(j)
        if max(MSroi2(:,i)) >= thresh
            imz=imz+1;
            mzroi_aug(imz)=mzroi2(1,i);
            MSroi_aug(:,imz)=[zeros(nr1,1);MSroi2(:,i)];
        else

            imzn=imzn+1;
            mzroi_non(imzn)=mzroi2(1,i);
            MSroi_non(:,imzn)=[zeros(nr1,1);MSroi2(:,i)];
        end
    end
end

[mzroi_aug_new,iorder] = sort(mzroi_aug);
mzroi_aug=mzroi_aug_new;
MSroi_aug=MSroi_aug(:,iorder);

if isfinite(mzroi_non)
    [mzroi_non_new,iorder] = sort(mzroi_non);
    mzroi_non=mzroi_non_new;
    MSroi_non=MSroi_non(:,iorder);
end