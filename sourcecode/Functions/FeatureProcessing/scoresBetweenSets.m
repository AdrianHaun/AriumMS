function [scoreArraySet1,scoreArraySet2] = scoresBetweenSets(set1,set2)
%% scoresBetweenSets calculates composite scores between each spectra of two sets
% Aligns both sets along the mass axis (rows), applies mass weighting 
% and calculates the composite score between each spectra. 
% 
% inputs: 
% set 1: matrix containing aligned spectra of set 1. column 1 must be the masses of
% each peak.
% set 2: matrix containing aligned spectra of set 2. column 1 must be the masses of
% each peak.
%
% outputs:
% scoreArraySet1: calculated composite score (column 1) and corresponding 
% spectra indices. Column 2: index of spectrum in set 1
% column 3: index of spectrum in set 2
%
% scoreArraySet2: calculated composite score (column 1) and corresponding 
% spectra indices. Column 2: index of spectrum in set 2
% column 3: index of spectrum in set 1

arguments
    set1 (:,:) double {mustBeNumeric,mustBeReal,mustHaveCorrectSize(set1)}
    set2 (:,:) double {mustBeNumeric,mustBeReal,mustHaveCorrectSize(set2)}
end

nSpectra_Set1 = width(set1)-1;
nSpectra_Set2 = width(set2)-1;

splitSpectraCells = cell(max([nSpectra_Set1;nSpectra_Set2]),2);

%split spectra into cells 
for iSpectrum = 1:max([nSpectra_Set1;nSpectra_Set2])
    if iSpectrum <= nSpectra_Set1
        splitSpectraCells{iSpectrum,1} = [set1(:,1),set1(:,iSpectrum+1)];
    end
    if iSpectrum <= nSpectra_Set2
        splitSpectraCells{iSpectrum,2} = [set2(:,1),set2(:,iSpectrum+1)];
    end
end
splitSpectraCells = reshape(splitSpectraCells,[],1);
splitSpectraCells(cellfun(@isempty,splitSpectraCells)) = [];

% align spectra
alignedSpectra = alignSpectra(splitSpectraCells,"normal","low","false");

%mass weighted intensities
alignedSpectra(:,2:end) = alignedSpectra(:,2:end).*alignedSpectra(:,1);
%remove mass column
alignedSpectra(:,1) = []; 
%re-scale
alignedSpectra = alignedSpectra./max(alignedSpectra);

%split into sets
set1 = alignedSpectra(:,1:nSpectra_Set1);
set2 = alignedSpectra(:,nSpectra_Set1+1:end);

%compare set1 to set2
scoreArraySet1 = zeros(nSpectra_Set1*nSpectra_Set2,3);
for iSpectrum_Set1 = 1:nSpectra_Set1
    for jSpectrum_Set2 = 1:nSpectra_Set2
        scoreArraySet1(iSpectrum_Set1*jSpectrum_Set2,1) = calculateCompositScore(set1(:,iSpectrum_Set1),set2(:,jSpectrum_Set2));
        scoreArraySet1(iSpectrum_Set1*jSpectrum_Set2,2:3) = [iSpectrum_Set1,jSpectrum_Set2];
    end
end

%compare set2 to set1
scoreArraySet2 = zeros(nSpectra_Set1*nSpectra_Set2,3);
for iSpectrum_Set2 = 1:nSpectra_Set2
    for jSpectrum_Set1 = 1:nSpectra_Set1
        scoreArraySet2(iSpectrum_Set2*jSpectrum_Set1,1) = calculateCompositScore(set2(:,iSpectrum_Set2),set1(:,jSpectrum_Set1));
        scoreArraySet2(iSpectrum_Set2*jSpectrum_Set1,2:3) = [iSpectrum_Set2,jSpectrum_Set1];
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
