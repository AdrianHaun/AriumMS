function msDataStruct = denoiseScans(msDataStruct,mode,cutOffValue)
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
    msDataStruct (1,1) struct
    mode        (1,1) string {mustBeMember(mode,["variable","threshold"])} = "variable";
    cutOffValue (1,1) double {mustBeFinite,mustBePositive} = 10
end

spectra = msDataStruct.spectra;

switch mode
    case "threshold"
        cutoff = cutOffValue;
    otherwise
        %% determine over all intensity bins
        ints = vertcat(spectra.rawScan);
        [~,edges] = histcounts(ints(:,2));
        cutoff = edges(2);
end

%% Clean Data
for j = 1:numel(spectra)
    data = spectra(j).rawScan;
    idx = data(:,2) <= cutoff;
    data(idx,2) = 0;
    spectra(j).processedScan = data;
end

msDataStruct.spectra = spectra;