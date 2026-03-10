function msDataStruct = compressScans(msDataStruct)
%% compressScans removes empty scans from MS-Data struct and compresses ScanData by removing entries without intensity
% DESCRIPTIVE TEXT

arguments
    msDataStruct (1,1) struct
end

spectra = msDataStruct.spectra;

%% Clean Data
for iSpectra = 1:numel(spectra)
    % compress processedScans
    processedData = spectra(iSpectra).processedScan;
    hasIntensity = processedData(:,2) > 0;
    spectra(iSpectra).processedScan = processedData(hasIntensity,:);
    % compress centroidScans
    centroidData = spectra(iSpectra).centroidedScan;
    hasIntensity = centroidData(:,2) > 0;
    spectra(iSpectra).centroidedScan = centroidData(hasIntensity,:);
end

% remove empty scans in spectra struct
% Determine which spectra have non-empty processedScan & centroidedScan
isEmptyScan = true(1,numel(spectra));
for iSpectra = 1:numel(spectra)
    % Consider a scan empty if rawScan is empty or missing, RT is NaN or msLevel is NaN
    hasProcessed = ~isempty(spectra(iSpectra).processedScan);
    hasCentroid = ~isempty(spectra(iSpectra).centroidedScan);
    if hasProcessed && hasCentroid
        isEmptyScan(iSpectra) = false;
    end
end
spectra = spectra(~isEmptyScan);

%update number of scans
nScans = numel(spectra);
msDataStruct.nSpectra       = nScans;

msDataStruct.spectra        = spectra;
