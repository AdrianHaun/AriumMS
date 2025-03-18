function [ScoresSet1,ScoresSet2] = ScoresBetweenSets(set1,set2)

arguments
    set1 (:,:) double {mustBeNumeric,mustBeReal,mustHaveCorrectSize(set1)}
    set2 (:,:) double {mustBeNumeric,mustBeReal,mustHaveCorrectSize(set2)}
end

nSpectra_Set1 = width(set1)-1;
nSpectra_Set2 = width(set2)-1;

splitSpectraCells = cell(max([nSpectra_Set1;nSpectra_Set2]),2);

%split spectra into cells
for n = 1:max([nSpectra_Set1;nSpectra_Set2])
    if n <= nSpectra_Set1
        splitSpectraCells{n,1} = [set1(:,1),set1(:,n+1)];
    end
    if n <= nSpectra_Set2
        splitSpectraCells{n,2} = [set2(:,1),set2(:,n+1)];
    end
end

splitSpectraCells = reshape(splitSpectraCells,[],1);
splitSpectraCells(cellfun(@isempty,splitSpectraCells)) = [];
% align spectra
alignedSpectra = AlignSpectra(splitSpectraCells,"normal");

%split into sets
alignedSpectra(:,1) = []; %remove mass column
set1 = alignedSpectra(:,1:nSpectra_Set1);
set2 = alignedSpectra(:,nSpectra_Set1+1:end);

%compare set1 to set2
ScoresSet1 = zeros(nSpectra_Set1*nSpectra_Set2,3);
for n = 1:nSpectra_Set1
    for m = 1:nSpectra_Set2
        ScoresSet1(n*m,1) = CompositScore(set1(:,n),set2(:,m));
        ScoresSet1(n*m,2:3) = [n,m];
    end
end

%compare set2 to set1
ScoresSet2 = zeros(nSpectra_Set1*nSpectra_Set2,3);
for n = 1:nSpectra_Set2
    for m = 1:nSpectra_Set1
        ScoresSet2(n*m,1) = CompositScore(set2(:,n),set1(:,m));
        ScoresSet2(n*m,2:3) = [n,m];
    end
end

end

% input validation function
function mustHaveCorrectSize(x)
% Test for more than 2 columns
if width(x) < 2
    eid = 'Size:incorrect';
    msg = 'Input must have at least two columns';
    error(eid,msg)
end
end