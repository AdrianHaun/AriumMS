function MsDataStruct = denoiseScans(MsDataStruct)
%DENOISESCANS Remove low-intensity peaks from mass spectra.
%
%   MsDataStruct = DENOISESCANS(MsDataStruct) determines an intensity
%   cutoff for each spectra in MsDataStruct.spectra and sets all
%   intensities below this cutoff to zero in the field processedScan.
%
%   Input:
%       MsDataStruct.spectra(i).rawScan - [m x 2] matrix with columns:
%           column 1: mass-to-charge (m/z)
%           column 2: intensity
%
%   Output:
%       MsDataStruct.spectra(i).processedScan - same size as rawScan,
%           with low-intensity values set to zero.
%

    arguments
        MsDataStruct (1,1) struct
    end

        % Input validation.
    assert(isfield(MsDataStruct, 'spectra'), ...
        'denoiseScans:MissingField', ...
        'Input struct must contain field ''spectra''.');

    assert(isfield(MsDataStruct, 'nSpectra'), ...
        'denoiseScans:MissingField', ...
        'Input struct must contain field ''nSpectra''.');

    spectra  = MsDataStruct.spectra;
    nSpectra = MsDataStruct.nSpectra;

    assert(isscalar(nSpectra) && isnumeric(nSpectra) && nSpectra >= 0 && ...
           nSpectra == floor(nSpectra), ...
        'denoiseScans:InvalidSpectraCount', ...
        'MsDataStruct.nSpectra must be a non-negative integer scalar.');

    assert(numel(spectra) == nSpectra, ...
        'denoiseScans:InconsistentSpectraCount', ...
        'MsDataStruct.nSpectra must equal numel(MsDataStruct.spectra).');

    for iSpectrum = 1:nSpectra
        assert(isfield(spectra(iSpectrum), 'rawScan'), ...
            'denoiseScans:MissingField', ...
            'spectra(%d) must contain field ''rawScan''.', iSpectrum);

        rawScan = spectra(iSpectrum).rawScan;

        % rawScan must be a numeric matrix with at least two columns.
        assert(isnumeric(rawScan) && ismatrix(rawScan) && size(rawScan, 2) >= 2, ...
            'denoiseScans:InvalidScanFormat', ...
            'spectra(%d).rawScan must be a numeric [m x 2] matrix [mass intensity ...].', ...
            iSpectrum);
    end

    % Determine overall intensity cutoff (max of second bin edge over spectra).
    cutoffPerSpectrum = zeros(nSpectra, 1);
    for iSpectrum = 1:nSpectra
        rawScan = spectra(iSpectrum).rawScan;
        intensities = full(rawScan(:, 2));
        [~, edges] = histcounts(intensities);
        cutoffPerSpectrum(iSpectrum) = edges(2);
    end


    % Apply cutoff to each spectrum.
    for iSpectrum = 1:nSpectra
        rawScan = full(spectra(iSpectrum).rawScan);
        lowIntensityMask = rawScan(:, 2) < cutoffPerSpectrum(iSpectrum);
        rawScan(lowIntensityMask, 2) = 0;
        spectra(iSpectrum).processedScan = rawScan;
    end

    MsDataStruct.spectra = spectra;
end