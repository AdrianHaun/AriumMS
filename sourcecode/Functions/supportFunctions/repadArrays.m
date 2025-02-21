function [ArrayOut,VecOut] = repadArrays(ArrayIn,VecIn,PadSize)
% Resizes array and vector with zeros
OrgSize=size(ArrayIn);
ArrayOut = padarray(ArrayIn,PadSize-OrgSize,0,'post');
VecOut = padarray(VecIn,PadSize(1)-OrgSize(1),0,'post');
end