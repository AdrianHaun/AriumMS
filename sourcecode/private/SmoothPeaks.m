function [MSroi,time] = SmoothPeaks(obj,MSroi,time)
Frame = obj.FrameSize;
Deg = obj.Degree;
parfor id = 1:size(MSroi,1)
    %depad Array
    oldSize=size(MSroi{id});
    MSroiTemp = MSroi{id};
    [MSroiTemp,timeTemp] = depadArrays(MSroiTemp,time{id});
    MSroiTemp = mssgolay(timeTemp,MSroiTemp,'Span',Frame,'Degree',Deg);
    %remove negativ, NaN and repad Array
    MSroiTemp=max(MSroiTemp,0);
    MSroiTemp(isnan(MSroiTemp))=0;
    [MSroi{id},time{id}] = repadArrays(MSroiTemp,timeTemp,oldSize);
    % set possible negative values to 0
    MSroi{id} = max(MSroi{id},0);
end

end
function [ArrayOut,VecOut] = depadArrays(ArrayIn,VecIn)
idx = VecIn~=0;
ArrayOut = ArrayIn(idx,:);
VecOut = VecIn(idx);
end

function [ArrayOut,VecOut] = repadArrays(ArrayIn,VecIn,PadSize)
OrgSize=size(ArrayIn);
ArrayOut = padarray(ArrayIn,PadSize-OrgSize,0,'post');
VecOut = padarray(VecIn,PadSize(1)-OrgSize(1),0,'post');
end