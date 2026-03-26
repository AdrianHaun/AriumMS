function MsDataStruct = compressScans(MsDataStruct)
% compressScans Remove empty scans from an MS data struct.
%   msDataStruct = compressScans(msDataStruct) removes peaks with
%   non-positive intensities from the processed and centroided scans of
%   each spectrum in msDataStruct.spectra. Spectra for which either
%   processedScan or centroidedScan is empty after compression are removed.
%
%   Input:
%       msDataStruct   Struct with fields:
%                        - spectra: 1-by-N struct array, each element
%                          containing at least fields processedScan and
%                          centroidedScan (double arrays with intensity in
%                          column 2).
%
%   Output:
%       msDataStruct   Same struct with compressed spectra and updated
%                      nSpectra field.

    arguments
        MsDataStruct (1,1) struct
    end

    spectraArray = MsDataStruct.spectra;

    %% Compress scan data by removing entries without intensity
    nSpectra = numel(spectraArray);

    for iSpectra = 1:nSpectra
        % Compress processedScan.
        processedData = spectraArray(iSpectra).processedScan;

        if ~isempty(processedData)
            if size(processedData, 2) < 2
                error('compressScans:InvalidProcessedScan', ...
                    ['processedScan for spectrum %d must have at least 2 ', ...
                     'columns (m/z and intensity).'], iSpectra);
            end

            hasIntensity = processedData(:, 2) > 0;
            spectraArray(iSpectra).processedScan = processedData(hasIntensity, :);
        end

        % Compress centroidedScan.
        centroidData = spectraArray(iSpectra).centroidedScan;

        if ~isempty(centroidData)
            if size(centroidData, 2) < 2
                error('compressScans:InvalidCentroidedScan', ...
                    ['centroidedScan for spectrum %d must have at least 2 ', ...
                     'columns (m/z and intensity).'], iSpectra);
            end

            hasIntensity = centroidData(:, 2) > 0;
            spectraArray(iSpectra).centroidedScan = centroidData(hasIntensity, :);
        end
    end

    %% Remove spectra with empty processed or centroided scans
    isEmptyScan = true(1, nSpectra);
    for iSpectra = 1:nSpectra
        hasProcessed = ~isempty(spectraArray(iSpectra).processedScan);
        hasCentroid = ~isempty(spectraArray(iSpectra).centroidedScan);

        % A scan is considered non-empty only if both are non-empty.
        if hasProcessed && hasCentroid
            isEmptyScan(iSpectra) = false;
        end
    end

    spectraArray = spectraArray(~isEmptyScan);

    % Update number of scans and spectra field.
    MsDataStruct.spectra  = spectraArray;
    MsDataStruct.nSpectra = numel(spectraArray);
end