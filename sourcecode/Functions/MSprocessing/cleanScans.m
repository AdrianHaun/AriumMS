function cleanedScans = cleanScans(rawScans)
%% cleanScans takes MS scans determines the lowest intensities and removes them
%
% Removes each intensity below the cutoff value. The cutoff value is 
% determined by histogram bin counts using the upper edge of bin one. 
%
% inputs: rawScans: cell array containing two column matrices, 
%                    column1: mass; column 2: intensity 
%
% output: cleanedScans: cleaned scans in the same format as the input

arguments
    rawScans (:,1) cell
end

cleanedScans = cell(size(rawScans));

%% determine over all intensity bins
ints = vertcat(rawScans{:});
[~,edges] = histcounts(ints(:,2));
cutoff = edges(2); 

%% Clean Data
parfor j = 1:height(cleanedScans)
    data = rawScans{j,1};
    idx = data(:,2) <= cutoff;
    data(idx,:) = [];
    cleanedScans{j,1} = data;
end
