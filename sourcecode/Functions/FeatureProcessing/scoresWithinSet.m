function scoreArray = scoresWithinSet(alignedSpectraMatrix)
%% scoresWithinSet calculates composite scores between each spectra of a set
% Applies mass weighting, then calculates the composite score between each 
% spectrum within the set 
% 
% inputs: 
% alignedSpectraMatrix: matrix containing aligned spectra. column 1 must be the masses of
% each peak.
%
% output:
% scoreArray: calculated composite score (column 1) and corresponding 
% spectra indices. Column 2: index of spectrum 1 column 3: index of spectrum 2
%

arguments
    alignedSpectraMatrix (:,:) double {mustBeNumeric,mustBeReal,mustHaveCorrectSize(alignedSpectraMatrix)}
end

intensities = alignedSpectraMatrix(:,2:end);
%mass weighted intensities
intensities = intensities.*alignedSpectraMatrix(:,1);

%scale each spectrum
intensities = intensities./max(intensities);

nSpectra = size(intensities,2);
if nSpectra > 1
    spectraIndex = 1:1:nSpectra;
    combinations = nchoosek(spectraIndex,2);
    scoreArray = zeros(size(combinations,1),1);
    for iCombination = 1:size(combinations,1)
        spectrum1 = intensities(:,combinations(iCombination,1))';
        spectrum2 = intensities(:,combinations(iCombination,2))';
        %remove entries that are zero in both spectra
        idx = spectrum1 == 0 & spectrum2 == 0;
        spectrum1(idx) = [];
        spectrum2(idx) = [];
        %calculate score
        scoreArray(iCombination,1) = calculateCompositeScore(spectrum1,spectrum2);
    end
    scoreArray = [scoreArray,combinations];
else
    scoreArray = [999,1,1];
end

% input validation function
function mustHaveCorrectSize(x)
% Test for more than 2 columns
if width(x) < 2
    eid = 'Size:incorrect';
    msg = 'Input must have at least two columns';
    error(eid,msg)
end