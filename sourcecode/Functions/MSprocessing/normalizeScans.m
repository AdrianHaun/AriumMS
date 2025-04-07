function normalizedScans = normalizeScans(originalScans)
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
    originalScans (:,1) cell
end

if isempty(originalScans)
    normalizedScans = cell(0,1);
    return
end

%preallocate output
normalizedScans = cell(size(originalScans));

parfor nScan = 1:height(normalizedScans)
    currentScan = originalScans{nScan,1};
    if isempty(currentScan)
        continue
    end
    currentScan(:,2) = currentScan(:,2)./max(currentScan(:,2));
    normalizedScans{nScan,1} = currentScan;
end