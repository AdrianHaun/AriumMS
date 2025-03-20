function [responses,deltaMZ] = ObjectiveFunc(output,targets,rtTol)

FeatIdentifier = output.FeatIdentifiers;
numFeat = size(FeatIdentifier,1);

if numFeat == 0
    responses = [0,1];
else
nTargets = size(targets,1);
mzTol = 0.05;
%check how many standard compounds are found
Found = zeros(height(targets),1);
deltaMZ = NaN(height(targets),1);
for n=1:height(targets)
    mzt=targets(n,1);
    RT=targets(n,2);
    idx=ismembertol(FeatIdentifier(:,1),mzt,mzTol,'DataScale',1);
    idxrt=ismembertol(FeatIdentifier(:,2),RT,rtTol,'DataScale',1);
    id=all([idx,idxrt],2);
    Found(n,1) = any(id);
    if sum(id)~=0
        deltaMZ(n,1) = min(abs(FeatIdentifier(id,1)-mzt));
    end
end
responses = [sum(Found)/nTargets,(numFeat-sum(Found))/numFeat];

end

responses = [responses*100,numFeat];
end