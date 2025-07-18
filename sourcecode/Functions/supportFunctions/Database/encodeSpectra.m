function encodedStrings = encodeSpectra(spectrumCells)
%% decodeSpectra decodes string mass spectra from database entries
% input: encodedSpectra as one dimensional string array
% output: cell array containing two column matrix. Column 1: mass, 
%           column 2: relative intensity

arguments
    spectrumCells (:,1) cell
end

encodedStrings = cell(size(spectrumCells,1),1);

parfor n = 1:height(encodedStrings)
    encoded = string(reshape(spectrumCells{n,1},1,[]));
    encodedStrings(n,1) = join(encoded,"-");
end