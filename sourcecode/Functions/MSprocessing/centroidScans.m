function centroidedScans = centroidScans(profileScans)
%% centroidScans takes profile MS scans and centroides them
%
%   Searches in each mass scan for mz peaks, takes all mz values within
%   half height and calculates the centroided mass as mean weighted by
%   intensity. Then sums all intensities used for mz calculation as
%   centroided intensity.
%
% inputs: profileScans: cell array containing two column matrices, 
%                    column1: mass; column 2: intensity 
%
% output: centroidedScans: centroided scans in the same format as the input

arguments
    profileScans (:,1) cell
end

if isempty(profileScans)
    centroidedScans = {[]};
    return
end

%preallocate output
centroidedScans = cell(size(profileScans));

parfor nScan = 1:height(centroidedScans)
    currentScan = profileScans{nScan,1};
    if isempty(currentScan)
        continue
    end
    %find mass peaks and width
    [maxIntensity,maxMZ,width] = findpeaks(currentScan(:,2),currentScan(:,1),'WidthReference','halfheight','SortStr','none');
    %preallocate centroided scan
    maxIntensity = zeros(height(maxIntensity),2);
    %compute weighted mean of mz and sum intensities
    for nPeak = 1:numel(width)
        id = abs(currentScan(:,1)-maxMZ(nPeak)) <= width(nPeak);
        vec = currentScan(id,:);
        mass = mean(vec(:,1),Weights = vec(:,2)/max(vec(:,2)));
        maxIntensity(nPeak,1) = mass; 
        maxIntensity(nPeak,2) = sum(vec(:,2));
    end
    centroidedScans{nScan,1} = maxIntensity;
end

%remove possible empty scans
idx = cellfun(@isempty,centroidedScans);
centroidedScans(idx) = [];