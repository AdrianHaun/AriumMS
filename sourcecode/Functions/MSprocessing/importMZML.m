function mzMLStruct = importMZML(fileLocation)
%IMPORTMZML Import mzML file and prepare structured mass-spectrometry data
%   mzMLStruct = IMPORTMZML(fileLocation) reads the mzML file specified by
%   fileLocation and returns a structured representation suitable for
%   downstream processing. The function performs the following high-level
%   steps (invoked after input validation):
%
%   1. read_mzml_mex - low-level mzML parsing that loads raw spectra and
%      metadata into an initial structure.
%   2. processMZMLStruct - extracts relevant data and structures it into
%      a final output struct.
%   3. denoiseScans - applies noise-reduction to MS1 and MS2 scans.
%   4. rmfield(...,"rawScan") - removes rawScan to free memory once denoised.
%   5. centroidScans - converts profile-mode spectra to centroided peaks.
%   6. normalizeMS2Scans - normalizes MS2 intensities across
%      spectra to facilitate comparisons.
%   7. compressScans - applies any compression/storage optimizations.
%
%   Input:
%       fileLocation - character vector specifying the full path to an
%                      existing .mzML file.
%
%   Output:
%       mzMLStruct - struct containing processed spectra, metadata, and
%                    auxiliary fields ready for analysis.
%
%   Notes:
%     - This function validates the input file using argument validation
%       (mustBeFile) and relies on auxiliary functions (read_mzml_mex,
%       processMZMLStruct, denoiseScans, centroidScans,
%       normalizeMS2Scans, compressScans) to perform the detailed work.
%     - Ensure the helper functions are on the MATLAB path before calling.
%
%   Example:
%       mz = importMZML("sample.mzML");
arguments (Input)
    fileLocation char {mustBeFile}
end

arguments (Output)
    mzMLStruct struct
end

mzMLStruct = read_mzml_mex(fileLocation);
mzMLStruct = processMZMLStruct(mzMLStruct);
mzMLStruct = denoiseScans(mzMLStruct);
mzMLStruct.spectra = rmfield(mzMLStruct.spectra,"rawScan");
mzMLStruct = centroidScans(mzMLStruct);
mzMLStruct = normalizeMS2Scans(mzMLStruct);
mzMLStruct = compressScans(mzMLStruct);

% throw error if file does not contain spectra data
if mzMLStruct.nSpectra == 0
    error('The mzML file does not contain any spectra data.');
end

% split MS1 and MS2
ms1Spectra = mzMLStruct.spectra([mzMLStruct.spectra.msLevel]' == 1);
ms2Spectra = mzMLStruct.spectra([mzMLStruct.spectra.msLevel]' == 2);

% Remove unnecessary fields from MS1 and MS2 spectra
ms1Spectra = rmfield(ms1Spectra, {'msLevel', 'precursorMz', 'collisionEnergy', 'fragmentationType', 'dataType'});
ms2Spectra = rmfield(ms2Spectra, {'msLevel', 'basePeakMz', 'basePeakIntensity', 'totalIonCurrent', 'dataType'});

% Remove old spectra field
mzMLStruct = rmfield(mzMLStruct,"spectra");

% Store split spectra
mzMLStruct.spectraMS1 = ms1Spectra;
mzMLStruct.spectraMS2 = ms2Spectra;