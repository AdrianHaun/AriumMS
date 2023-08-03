function CorrectedPeakData = CWTBorderCorrection(PeakData,EIC,smoothedEIC)
%friction border correction
smoothedEIC=rescale(smoothedEIC,0,1);
IntThresh=0.001;
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
        if EIC(u-1)-mean(EIC(u-1:u+7)) <= std(EIC(u-1:u+7)) || round(EIC(u-1)-EIC(u))==0
            u=u-1;
        else
            PeakData(n,3)=u;
            break
        end
    end
    l=PeakData(n,2);
    while l<=numel(EIC)-1 && l>8 %lower bond
        if EIC(l+1)-mean(EIC(l-7:l+1)) <= std(EIC(l-7:l+1)) || round(EIC(l)-EIC(l+1))==0
            l=l+1;
        else
            PeakData(n,2)=l;
            break
        end
    end
end
CorrectedPeakData = PeakData;
end