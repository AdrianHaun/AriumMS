function MzMLStruct = importMZML(fileLocation)
%IMPORTMZML Import mzML file and prepare structured mass-spectrometry data.
%   MzMLStruct = IMPORTMZML(fileLocation) reads the mzML file specified by
%   fileLocation and returns a structured representation suitable for
%   downstream processing. The function performs the following high-level
%   steps (invoked after input validation):
%
%   1. read_mzml_mex      - low-level mzML parsing that loads raw spectra
%                           and metadata into an initial structure.
%   2. processMZMLStruct  - extracts relevant data and structures it into
%                           a final output struct.
%   3. denoiseScans       - applies noise reduction to MS1 and MS2 scans.
%   4. rmfield(...,'rawScan')
%                        - removes rawScan to free memory once denoised.
%   5. centroidScans      - converts profile-mode spectra to centroided
%                           peaks.
%   6. normalizeMS2Scans  - normalizes MS2 intensities across spectra to
%                           facilitate comparisons.
%   7. compressScans      - applies any compression/storage optimizations.
%
%   Input:
%       fileLocation - character vector specifying the full path to an
%                      existing .mzML file.
%
%   Output:
%       MzMLStruct  - struct containing processed spectra, metadata, and
%                     auxiliary fields ready for analysis.
%
%   Notes:
%     - This function validates the input file using argument validation
%       (mustBeFile) and relies on auxiliary functions (read_mzml_mex,
%       processMZMLStruct, denoiseScans, centroidScans,
%       normalizeMS2Scans, compressScans) to perform the detailed work.
%     - Ensure the helper functions are on the MATLAB path before calling.
%
%   Example:
%       mz = importMZML('sample.mzML');

arguments (Input)
    fileLocation char {mustBeFile}
end

arguments (Output)
    MzMLStruct struct
end

% Load and process raw mzML data.
MzMLStruct = read_mzml_mex(fileLocation);
MzMLStruct = processMZMLStruct(MzMLStruct);
MzMLStruct = denoiseScans(MzMLStruct);

% Remove large raw scan data to free memory.
if isfield(MzMLStruct, 'spectra') && isfield(MzMLStruct.spectra, 'rawScan')
    MzMLStruct.spectra = rmfield(MzMLStruct.spectra, 'rawScan');
end

MzMLStruct = centroidScans(MzMLStruct);
MzMLStruct = normalizeMS2Scans(MzMLStruct);
MzMLStruct = compressScans(MzMLStruct);

% Throw error if file does not contain spectra data.
if ~isfield(MzMLStruct, 'nSpectra') || MzMLStruct.nSpectra == 0
    error('importMZML:NoSpectra', ...
        'The mzML file does not contain any spectra data.');
end

% Split MS1 and MS2 spectra.
if ~isfield(MzMLStruct, 'spectra') || isempty(MzMLStruct.spectra)
    error('importMZML:MissingSpectraField', ...
        'MzMLStruct.spectra is missing or empty after processing.');
end

msLevels   = [MzMLStruct.spectra.msLevel]';
ms1Mask    = (msLevels == 1);
ms2Mask    = (msLevels == 2);

ms1Spectra = MzMLStruct.spectra(ms1Mask);
ms2Spectra = MzMLStruct.spectra(ms2Mask);

% Remove unnecessary fields from MS1 and MS2 spectra.
if ~isempty(ms1Spectra)
    ms1Spectra = rmfield(ms1Spectra, ...
        {'msLevel','precursorMz','collisionEnergy', ...
         'fragmentationType','dataType'});
end

if ~isempty(ms2Spectra)
    ms2Spectra = rmfield(ms2Spectra, ...
        {'msLevel','basePeakMz','basePeakIntensity', ...
         'totalIonCurrent','dataType'});
end

% Remove old spectra field and store split spectra.
MzMLStruct = rmfield(MzMLStruct, 'spectra');
MzMLStruct.spectraMS1 = ms1Spectra;
MzMLStruct.spectraMS2 = ms2Spectra;

end