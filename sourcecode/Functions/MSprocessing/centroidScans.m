function msDataStruct = centroidScans(msDataStruct)
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
    msDataStruct (1,1) struct
end

if isempty(msDataStruct)
    return
end

spectras = msDataStruct.spectra;

parfor iScan = 1:msDataStruct.nSpectra

    %check if already centroid
    if spectras(iScan).dataType == "profile"
        currentScan = spectras(iScan).processedScan;

        %find mass peaks and width
        [maxIntensity,maxMZ,width] = findpeaks(currentScan(:,2),currentScan(:,1),'WidthReference','halfheight','SortStr','none');
        %preallocate centroided scan
        maxIntensity = zeros(height(maxIntensity),2);
        %compute weighted mean of mz and sum intensities
        for jPeak = 1:numel(width)
            id = abs(currentScan(:,1)-maxMZ(jPeak)) <= width(jPeak);
            vec = currentScan(id,:);
            mass = mean(vec(:,1),Weights = vec(:,2)/max(vec(:,2)));
            maxIntensity(jPeak,1) = mass;
            maxIntensity(jPeak,2) = max(vec(:,2));
        end
        spectras(iScan).centroidedScan = maxIntensity;

    else
        spectras(iScan).centroidedScan = spectras(iScan).processedScan;
    end
end

%store
msDataStruct.spectra = spectras;