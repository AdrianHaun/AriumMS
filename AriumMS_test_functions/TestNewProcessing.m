function [responses,deltaMZ] = TestNewProcessing(output,rtTol)

targetDataTable = load("TargetDataNew.mat");
targetDataTable = targetDataTable.targetDataTable;

targetCells = cell(size(output));
for n = 1:length(output)
%build target value array
name = output(n).groupName;
vals = [targetDataTable.M,targetDataTable.(name)];
vals(vals(:,2)==0,:) = [];
targetCells{n} = vals;
end


responses = zeros(length(targetCells),3);

for f = 1:length(targetCells)

targets = targetCells{f};

FeatIdentifier = [output(f).feature.mass_measured;output(f).feature.retentionTime]';
numFeat = size(FeatIdentifier,1);

if numFeat == 0
    out = [0,1];
else

    nTargets = size(targets,1);
    mzTol = 0.01;
    %check how many standard compounds are found
    Found = zeros(height(targets),1);
    deltaMZ = NaN(height(targets),1);
    for n = 1:height(targets)
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
    out = [sum(Found)/nTargets,(numFeat-sum(Found))/numFeat];

end

responses(f,:) = [out*100,numFeat];
end