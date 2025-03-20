function [responses,deltaMZ] = ObjectiveTest(output,targets,rtTol)

FeatIdentifier = output.FeatIdentifiers;
numFeat = size(FeatIdentifier,1);

if numFeat == 0
    responses = [0,1];
else
nTargets = size(targets,1);
mzTol = 0.03;
%check how many standard compounds are found
Found = zeros(height(targets),1);
deltaMZ = NaN(height(targets),1);
for n=1:height(targets)
    mzt=targets(n,1);
    RT=targets(n,2);
    idx = abs(FeatIdentifier(:,1) - mzt)<=mzTol;
    idxrt = abs(FeatIdentifier(:,2) - RT)<=rtTol;
    id=all([idx,idxrt],2);
    Found(n,1) = any(id);
    if sum(id)~=0
        deltaMZ(n,1) = (min(abs(FeatIdentifier(id,1)-mzt))/mzt*1e6);
    end
end
responses = [sum(Found)/nTargets,(numFeat-sum(Found))/numFeat];

end

responses = [responses*100,numFeat];
end