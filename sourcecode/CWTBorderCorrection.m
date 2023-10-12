function CorrectedPeakData = CWTBorderCorrection(PeakData,EIC,smoothedEIC)
% Performs friction Border correction and Moving Standard Deviation Border
% Correction on initial Peak borders from AutoCWT

EIC=full(EIC);
smoothedEIC = full(smoothedEIC);
noise=EIC-smoothedEIC;
%friction border correction
smoothedEIC=rescale(smoothedEIC,0,1);
IntThresh=0.0001;
%normalize thresh
IntThresh=(max(smoothedEIC)-min(smoothedEIC))*IntThresh;
for n=1:size(PeakData,1)
    u=PeakData(n,3);
    while u<numel(smoothedEIC)-1 %upper bond
        if smoothedEIC(u)-smoothedEIC(u+1)>IntThresh
            u=u+1;
        else
            PeakData(n,3)=u;
            break
        end
    end
    l=PeakData(n,2);
    while l>1 %lower bond
        if smoothedEIC(l)-smoothedEIC(l-1)>IntThresh
            l=l-1;
        else
            PeakData(n,2)=l;
            break
        end
    end
end
%moving STD correction to remove errors from smoothed peaks
for n=1:size(PeakData,1)
    u=PeakData(n,3);
    while u>1 && u<numel(EIC)-8 %upper bond
        if abs(EIC(u-1)-mean(EIC(u-1:u+7,1))) <= std(noise(u-1:u+7,1)) || round(EIC(u-1,1)-EIC(u,1))==0
            u=u-1;
        else
            PeakData(n,3)=u;
            break
        end
    end
    l=PeakData(n,2);
    while l<=numel(EIC)-1 && l>8 %lower bond
        if abs(EIC(l+1,1)-mean(EIC(l-7:l+1,1))) <= std(noise(l-7:l+1,1)) || round(EIC(l,1)-EIC(l+1,1))==0
            l=l+1;
        else
            PeakData(n,2)=l;
            break
        end
    end
end
%get final peak location and height
for n=1:size(PeakData,1)
    EICtemp = EIC;
    EICtemp(1:PeakData(n,2)-1)=0;
    EICtemp(PeakData(n,3)+1:end)=0;
    [PeakData(n,4),PeakData(n,1)]=max(EICtemp);
end
CorrectedPeakData = PeakData;
end