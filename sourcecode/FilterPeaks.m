function [peaks,tempStorage] = FilterPeaks(peaks,MinPWDataPoints,MaxPWDataPoints,maxSN,Noise,EIC)
% Filters identified peaks from AutoCWT
                tempStorage = cell(8,1);
                %% Peak filter
                %remove peaks with wrong boundaries
                if isempty(peaks)==false
                    idx=peaks(:,2)>=peaks(:,3);
                    peaks(idx,:)=[];
                end
                %remove duplicate peaks
                if isempty(peaks)==false
                    peaks=unique(peaks,'rows');
                end
                %remove peaks with height = 0
                if isempty(peaks)==false
                    idx=peaks(:,4)==0;
                    peaks(idx,:)=[];
                end
                %remove peaks with bad Peak asymmetry
                if isempty(peaks)==false
                    idx= (peaks(:,3)-peaks(:,1))./(peaks(:,1)-peaks(:,2));
                    idx = idx<0.3 | idx>3;
                    peaks(idx,:)=[];
                end
                %less than minimum peak width
                if isempty(peaks)==false
                    idx=peaks(:,3)-peaks(:,2)<MinPWDataPoints;
                    peaks(idx,:)=[];
                    tempStorage{5,1}=sum(idx);
                end
                %more than maximum peak width
                if isempty(peaks)==false
                    idx=peaks(:,3)-peaks(:,2)>MaxPWDataPoints;
                    peaks(idx,:)=[];
                    tempStorage{6,1}=sum(idx);
                end
                %S/N peak rejection
                SN=zeros(size(peaks,1),1);
                if isempty(peaks)==false
                    SN=peaks(:,4)./Noise;
                    idx=SN<maxSN;
                    peaks(idx,:)=[];
                    SN(idx,:)=[];
                    tempStorage{7,1}=sum(idx);
                end
                % entropy calculation
                 if isempty(peaks)==false
                    Entropy = CalculatePeakEntropy(peaks,full(EIC));
                 else
                    Entropy = zeros(length(SN),1);
                 end
                tempStorage{4,1}=[Entropy,SN];
    end