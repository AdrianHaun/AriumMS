function [ArrayOut,VecOut] = depadArrays(ArrayIn,VecIn)
% Removes padding from Arrays
idx = VecIn~=0;
ArrayOut = ArrayIn(idx,:);
VecOut = VecIn(idx);
end