function denoisedScans = denoiseScans(noisyScans)
%% denoiseScans takes MS scans determines the noise level and sets noise to zero
%
% Sets each intensity below the cutoff value to zero. The cutoff value is 
% determined by histogram bin counts using the upper edge of bin one. 
%
% inputs: profileScans: cell array containing two column matrices, 
%                    column1: mass; column 2: intensity 
%
% output: denoisedScans: denoised scans in the same format as the input

arguments
    noisyScans (:,1) cell
end

%preallocate output
denoisedScans = cell(size(noisyScans));

%calculate overall noise level
ints = vertcat(noisyScans{:});
[~,edges] = histcounts(ints(:,2));
cutoff = edges(2); 

parfor nScan = 1:height(denoisedScans)
    currentScan = noisyScans{nScan,1};
    if isempty(currentScan)
        continue
    end
    currentScan(currentScan(:,2) <= cutoff,2) = 0;
    denoisedScans{nScan,1} = currentScan;
end