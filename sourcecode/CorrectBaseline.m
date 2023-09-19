function [MSroi,time] = CorrectBaseline(obj,MSroi,time)
WSize = obj.WindowSize;
SSize = obj.StepSize;
RegMethod = obj.RegressionMethod;
EstMethod =obj.EstimationMethod;
SmooMethod = obj.SmoothMethod;
Quan = obj.QuantilVal;
parfor id = 1:size(MSroi,1)
    oldSize=size(MSroi{id});
    %depad Array
    MSroiTemp = MSroi{id};
    [MSroiTemp,timeTemp] = depadArrays(MSroiTemp,time{id});
    MSroiTemp = msbackadj(timeTemp,MSroiTemp,'WindowSize',WSize,'StepSize',SSize,'RegressionMethod',RegMethod,'EstimationMethod',EstMethod,'SmoothMethod',SmooMethod,'QuantileValue',Quan,'PreserveHeights',true);
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