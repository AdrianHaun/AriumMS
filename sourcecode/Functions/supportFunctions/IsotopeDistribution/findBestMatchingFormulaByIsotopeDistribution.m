function bestIndex = findBestMatchingFormulaByIsotopeDistribution(decomposition,measuredIsotopePattern)

% convert decomposition to structure


decomposition = cell2struct(num2cell(double(decomposition), 1), {'C' 'H' 'Br' 'Cl' 'F' 'I' 'N' 'O' 'P' 'S'}, 2);

% calculate theoretical distros
[calculatedDistributions, ~, ~] = isotopicdist(decomposition,"NoiseThreshold",100);

%compare
for iDecomposition = 1:height(calculatedDistributions)
    
    calculatedDistro = calculatedDistributions{iDecomposition};
    %align and normalize
    alignedDistros = alignSpectra(spectraCells,"normal","medium","true",0.025);
    %compare
    
end