function PeakEntropy = CalculatePeakEntropy(PeakData,EIC)
% Calculates Peak entropy for all peaks
D=diff(EIC);
p=zeros(size(PeakData.peakLocation));
for n=1:numel(p)
    %extract peak range
    Peak = D(PeakData.peakStartLocation(n,:):PeakData.peakEndLocation(n,:));
    maxidx = PeakData.peakLocation(n)-PeakData.peakStartLocation(n);
    % check normal or variant point , variant point = 1
    premax=Peak(1:maxidx-1)<0;
    postmax=Peak(maxidx+1:end)>0;
    VarPoints=[premax; false; postmax];
    p(n,1)=sum(VarPoints)/numel(VarPoints);
end
PeakEntropy=-p.*log2(p)-(1-p).*log2(1-p);
PeakEntropy(isnan(PeakEntropy))=0;
end