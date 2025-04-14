function cleanedScans = denoiseScans(rawScans,mode,cutOffValue)
%% cleanScans takes MS scans determines the lowest intensities and removes them
%
% Removes each intensity below the cutoff value. The cutoff value is
% determined by histogram bin counts using the upper edge of bin one.
%
% inputs: rawScans: cell array containing two column matrices,
%                    column1: mass; column 2: intensity
%         mode: switches between variable mode (noise level is determined by histcounts)
%               or threshold mode (cutoff value is set by the user)
%         cutOffValue: user defined intensity cutoff when using threshold mode
%
% output: cleanedScans: cleaned scans in the same format as the input

arguments
    rawScans    (:,1) cell
    mode        (1,1) string {mustBeMember(mode,["variable","threshold"])} = "variable";
    cutOffValue (1,1) double {mustBeFinite,mustBePositive} = 0.01
end

if isempty(rawScans)
    cleanedScans = cell(0,1);
    return
end

cleanedScans = cell(size(rawScans));

switch mode
    case "threshold"
        cutoff = cutOffValue;
    otherwise
        %% determine over all intensity bins
        ints = vertcat(rawScans{:});
        [~,edges] = histcounts(ints(:,2));
        cutoff = edges(2);
end

%% Clean Data
parfor j = 1:height(cleanedScans)
    data = rawScans{j,1};
    idx = data(:,2) <= cutoff;
    data(idx,:) = [];
    cleanedScans{j,1} = data;
end
