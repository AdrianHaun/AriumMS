function msDataStruct = normalizeMS2Scans(msDataStruct)
%% normalizeScans takes MS scans and normalizes each scan 
%
% The intensities within each scan are normalized to to maximum intensity
% in each scan. Output intensities therefor range from 0 to 1.
%
% inputs: originalScans: cell array containing two column matrices, 
%                    column1: mass; column 2: intensity 
%
% output: normalizedScans: normalized scans in the same format as the input

arguments
    msDataStruct (1,1) struct
end

spectras = msDataStruct.spectra;

for nScan = 1:numel(spectras)
    currentScan = spectras(nScan).centroidedScan;
    if spectras(nScan).msLevel == 2
        currentScan(:,2) = currentScan(:,2)./max(currentScan(:,2));
        spectras(nScan).centroidedScan = currentScan;
    else
        continue
    end
end

msDataStruct.spectra = spectras;