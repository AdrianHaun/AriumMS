function decodedSpectrum = decodeSpectra(encodedSpectrum)
%% decodeSpectra decodes string mass spectra from database entries
% input: encodedSpectra as one dimensional string array
% output: cell array containing two column matrix. Column 1: mass, 
%           column 2: relative intensity

arguments
    encodedSpectrum (:,1) string
end

decodedSpectrum = cell(size(encodedSpectrum,1),1);

parfor n = 1:size(encodedSpectrum,1)
    decoded = split(encodedSpectrum(n,1),"-",2);
    decodedSpectrum{n,1} = str2double(reshape(decoded,[],2));
end