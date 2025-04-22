function NumComponents = DetermineNumComponents(Spectras)
%SVD
SingularVals = svd(Spectras);
%Determine number of independent vectors that explain 95% of variance
SingularVals = SingularVals./sum(SingularVals);
SingularVals = cumsum(SingularVals);
aRef = 0.95; %Set reference
aDiff = abs(SingularVals - aRef); %Calculate diff
[~, NumComponents] = min(aDiff); %Find value closest to aRef
end