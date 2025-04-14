function decodedSpectrum = decodeSpectra(encodedSpectrum)
%% decodeSpectra decodes base64 encoded mass spectra

% input: encodedSpectra as one dimensional string array
% output: cell array containing two column matrix. Column 1: mass, 
%           column 2: relative intensity

arguments
    encodedSpectrum (:,1) string
end

decodedSpectrum = cell(size(encodedSpectrum,1),1);

parfor n = 1:size(encodedSpectrum,1)
    decoded = matlab.net.base64decode(encodedSpectrum(n,1));
    decoded = typecast(decoded,'double');
    decodedSpectrum{n,1} = reshape(decoded,[],2);
end