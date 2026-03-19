function msDataStruct = denoiseScans(msDataStruct)
%% cleanScans takes MS scans determines the lowest intensities and removes them
%
% Removes each intensity below the cutoff value. The cutoff value is
% determined by histogram bin counts using the upper edge of bin one.
%
% inputs: rawScans: cell array containing two column matrices,
%                    column1: mass; column 2: intensity
%


arguments
    msDataStruct (1,1) struct
end

spectra = msDataStruct.spectra;

nSpectra = msDataStruct.nSpectra;

%% determine over all intensity bins
cutoff = zeros(nSpectra,1);
for iSpectrum = 1:nSpectra
    [~,edges] = histcounts(full(spectra(iSpectrum).rawScan(:,2)));
    cutoff(iSpectrum,1) = edges(2);
end
cutoff = max(cutoff);

%% Clean Data
for iSpectrum = 1:nSpectra
    data = full(spectra(iSpectrum).rawScan);
    idx = data(:,2) < cutoff;
    data(idx,2) = 0;
    spectra(iSpectrum).processedScan = data;
end

msDataStruct.spectra = spectra;