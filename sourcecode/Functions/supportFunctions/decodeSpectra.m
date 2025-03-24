function decodedSpectrum = decodeSpectra(encodedSpectra)
%% decodeSpectra decodes base64 encoded mass spectra

% input: encodedSpectra as one dimensional string array
% output: cell array containing two column matrix. Column 1: mass, 
%           column 2: relative intensity

arguments
    encodedSpectra (:,1) string
end

decodedSpectrum = cell(size(encodedSpectra,1),1);

parfor n = 1:size(encodedSpectra,1)
    decoded = matlab.net.base64decode(encodedSpectra(n,1));
    decoded = typecast(decoded,'double');
    decodedSpectrum{n,1} = reshape(decoded,[],2);
end