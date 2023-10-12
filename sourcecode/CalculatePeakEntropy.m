function PeakEntropy = CalculatePeakEntropy(PeakData,EIC)
% Calculates Peak entropy for all peaks
D=diff(EIC);
p=zeros(size(PeakData,1),1);
for n=1:size(PeakData,1)
    %extract peak range
    Peak = D(PeakData(n,2):PeakData(n,3));
    maxidx=PeakData(n,1)-PeakData(n,2);
    % check normal or variant point , variant point = 1
    premax=Peak(1:maxidx-1)<0;
    postmax=Peak(maxidx+1:end)>0;
    VarPoints=[premax; false; postmax];
    p(n,1)=sum(VarPoints)/numel(VarPoints);
end
PeakEntropy=-p.*log2(p)-(1-p).*log2(1-p);
PeakEntropy(isnan(PeakEntropy))=0;
end