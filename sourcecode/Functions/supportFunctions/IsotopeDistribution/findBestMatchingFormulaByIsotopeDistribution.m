function scores = findBestMatchingFormulaByIsotopeDistribution(decomposition,measuredIsotopePattern,resolutionDa)

% convert decomposition to structure
decomposition = cell2struct(num2cell(double(decomposition), 1), {'C' 'H' 'Br' 'Cl' 'F' 'I' 'N' 'O' 'P' 'S'}, 2);

% calculate theoretical distros
[calculatedDistributions, ~, ~] = isotopicdist(decomposition,"NoiseThreshold",100,"Resolution",resolutionDa);
% pack into cell if only one decomposition
if ~iscell(calculatedDistributions)
    calculatedDistributions = {calculatedDistributions};
end

scores = zeros(height(decomposition),1);
%compare
for iDecomposition = 1:height(calculatedDistributions)
    calculatedDistro = calculatedDistributions{iDecomposition};
    %align and normalize
    alignedDistros = alignSpectra([{calculatedDistro};{measuredIsotopePattern}],"normal","medium","true",0.025);
    %compare
    scores(iDecomposition,1) = calculateCompositeScore(alignedDistros(:,2),alignedDistros(:,3));
end