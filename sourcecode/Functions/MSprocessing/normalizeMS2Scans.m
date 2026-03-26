function MsDataStruct = normalizeMS2Scans(MsDataStruct)
% normalizeMS2Scans Normalizes intensities of MS2 scans in an MS data struct.
%
% The intensities within each MS2 scan (msLevel == 2) are normalized to the
% maximum intensity in that scan. Output intensities therefore range from 0 to 1.
%
% Input:
%   msDataStruct  Struct with field 'spectra'. Each element of spectra is
%                 expected to contain:
%                   - msLevel        : numeric MS level
%                   - centroidedScan : two-column matrix [mass, intensity]
%
% Output:
%   msDataStruct  Same struct as input, but with centroidedScan(:,2) for
%                 MS2 scans normalized to their per-scan maximum.
%
%   msDataStruct = normalizeMS2Scans(msDataStruct)

arguments
    MsDataStruct (1,1) struct
end

spectraArray = MsDataStruct.spectra;

for iScan = 1:numel(spectraArray)
    if spectraArray(iScan).msLevel == 2
        currentScan = spectraArray(iScan).centroidedScan;

        % Normalize intensity column by its maximum value.
        currentScan(:, 2) = currentScan(:, 2) ./ max(currentScan(:, 2));

        spectraArray(iScan).centroidedScan = currentScan;
    end
end

MsDataStruct.spectra = spectraArray;
end