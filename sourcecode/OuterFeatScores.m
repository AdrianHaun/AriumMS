function BetweenGroupScores = OuterFeatScores(MSnSpectraCells,mzTol,TolUnit)

CanCompare = sum(~cellfun("isempty",MSnSpectraCells),2)>=2;
MSnSpectraCells = cellfun(@(x) double(x),MSnSpectraCells,'UniformOutput',false);
NumCombis = size(nchoosek(1:1:size(MSnSpectraCells,2),2),1);
BetweenGroupScores = cell(size(MSnSpectraCells,1),NumCombis);
for n = 1:size(MSnSpectraCells,1)
    if CanCompare(n) == false
        continue
    end
    Spectras = MSnSpectraCells(n,:);
    nSpectras = cellfun(@(x) size(x,2)-1,Spectras);

    %Merge into Matrix
    Spectras = mergeMatricesWithTolerance(Spectras',mzTol,TolUnit);

    Spectras(:,1) = [];
    %scale each spectrum
    Spectras = Spectras./max(Spectras);
    %split matrix into cell array
    Spectras = mat2cell(Spectras,size(Spectras,1),nSpectras);
    % Get the number of matrices in the cell array
    numGroups = numel(Spectras);
    GroupCombinations = nchoosek(1:1:numGroups,2);
    GroupScores = cell(1,size(GroupCombinations,1));
    for i = 1:size(GroupCombinations,1)
        Group1 = Spectras{GroupCombinations(i,1)};
        Group2 = Spectras{GroupCombinations(i,2)};
        nSpectra1 = size(Group1,2);
        nSpectra2 = size(Group2,2);
        SpectraCombinations = table2array(combinations(1:1:nSpectra1,1:1:nSpectra2));
        ScoreVec = zeros(size(SpectraCombinations,1),1);
        for j = 1:size(SpectraCombinations,1)
            
            Spectra1 = Group1(:,SpectraCombinations(j,1));
            Spectra2 = Group2(:,SpectraCombinations(j,2));
            idx = Spectra1 == 0 & Spectra2 == 0;
            Spectra1(idx) = [];
            Spectra2(idx) = [];
            ScoreVec(j) = IdentityMatchFactor(Spectra1', Spectra2');
        end
        ScoreVec = [ScoreVec,SpectraCombinations];
        GroupScores{i} = ScoreVec;
    end
    BetweenGroupScores(n,:) = GroupScores;
end