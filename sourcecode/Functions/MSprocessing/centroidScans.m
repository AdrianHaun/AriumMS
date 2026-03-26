function MsDataStruct = centroidScans(MsDataStruct)
%CENTROIDSCANS Centroid profile mass spectra.
%
%   MsDataStruct = CENTROIDSCANS(MsDataStruct) takes profile-mode mass
%   spectra stored in MsDataStruct.spectra(i).processedScan and computes
%   centroided spectra. For each profile spectrum, peaks are found using
%   FINDPEAKS on intensity versus m/z. For each peak, all points within
%   the half-height width are collected, the centroid m/z is computed as
%   an intensity-weighted mean, and the centroid intensity is the maximum
%   intensity in that window. Non-profile spectra are passed through
%   unchanged.
%
%   Input:
%       MsDataStruct.spectra(i).processedScan - [m x 2] matrix:
%           column 1: m/z
%           column 2: intensity
%       MsDataStruct.spectra(i).dataType      - "profile" or other.
%
%   Output:
%       MsDataStruct.spectra(i).centroidedScan - [k x 2] matrix:
%           column 1: centroid m/z
%           column 2: centroid intensity
%

arguments
    MsDataStruct (1,1) struct
end

% Early exit for empty struct.
if isempty(MsDataStruct)
    return
end

%% Input Validation
assert(isfield(MsDataStruct, 'spectra'), ...
    'centroidScans:MissingField', ...
    'Input struct must contain field ''spectra''.');

assert(isfield(MsDataStruct, 'nSpectra'), ...
    'centroidScans:MissingField', ...
    'Input struct must contain field ''nSpectra''.');

spectra  = MsDataStruct.spectra;
nSpectra = MsDataStruct.nSpectra;

assert(isscalar(nSpectra) && isnumeric(nSpectra) && nSpectra >= 0 && ...
    nSpectra == floor(nSpectra), ...
    'centroidScans:InvalidSpectraCount', ...
    'MsDataStruct.nSpectra must be a non-negative integer scalar.');

assert(numel(spectra) == nSpectra, ...
    'centroidScans:InconsistentSpectraCount', ...
    'MsDataStruct.nSpectra must equal numel(MsDataStruct.spectra).');

for iScan = 1:nSpectra
    assert(isfield(spectra(iScan), 'processedScan'), ...
        'centroidScans:MissingField', ...
        'spectra(%d) must contain field ''processedScan''.', iScan);

    scan = spectra(iScan).processedScan;
    assert(isnumeric(scan) && ismatrix(scan) && size(scan, 2) >= 2, ...
        'centroidScans:InvalidScanFormat', ...
        'spectra(%d).processedScan must be a numeric [m x 2+] matrix [m/z intensity ...].', ...
        iScan);

    assert(isfield(spectra(iScan), 'dataType'), ...
        'centroidScans:MissingField', ...
        'spectra(%d) must contain field ''dataType''.', iScan);
end

%% Main processing loop.
for iScan = 1:nSpectra

    % check if input scan is profile
    if strcmp(spectra(iScan).dataType, 'profile')

        currentScan = spectra(iScan).processedScan;
        mz  = currentScan(:, 1);
        int = currentScan(:, 2);

        % Find peaks and half-height widths.
        [peakIntensity, peakMZ, width] = findpeaks( ...
            int, mz, ...
            'WidthReference', 'halfheight', ...
            'SortStr',        'none');

        nPeaks = numel(peakIntensity);
        centroidedScan = zeros(nPeaks, 2);  % [m/z, intensity]

        % For each peak, compute centroid.
        for jPeak = 1:nPeaks
            % Define window using m/z bounds (assumes mz sorted).
            leftBound  = peakMZ(jPeak) - width(jPeak);
            rightBound = peakMZ(jPeak) + width(jPeak);

            inWindow = (mz >= leftBound) & (mz <= rightBound);

            mzWindow  = mz(inWindow);
            intWindow = int(inWindow);

            if isempty(mzWindow)
                % Should not really happen, but guard for safety.
                centroidedScan(jPeak, 1) = peakMZ(jPeak);
                centroidedScan(jPeak, 2) = peakIntensity(jPeak);
            else
                % Intensity-weighted centroid: sum(mz*I)/sum(I).
                sumWeights     = sum(intWindow);
                weightedMzSum  = sum(mzWindow .* intWindow);
                massCentroid   = weightedMzSum / sumWeights;

                localMax       = max(intWindow);

                centroidedScan(jPeak, 1) = massCentroid;
                centroidedScan(jPeak, 2) = localMax;
            end
        end

        spectra(iScan).centroidedScan = centroidedScan;

    else
        spectra(iScan).centroidedScan = spectra(iScan).processedScan;
    end
end

MsDataStruct.spectra = spectra;
end