function peakXIC = extractPeakXIC(xic,peakBorder)
%% extractPeakXIC gets XIC of peaks from cell array of file XIC
% Inputs:   xic: cell array of extracted ion chromatograms
%           peakBorder: cell array of peak borders

arguments
    xic (:,1) cell
    peakBorder (:,:) cell
end

%reshape input to column vector
[originalnRows,originalnCols] = size(peakBorder);
peakBorder = reshape(peakBorder,[],1);
%preallocation
peakXIC = cell(size(peakBorder));
for iPeak = 1:size(peakXIC,1)
    lowerBorder = peakBorder{iPeak,1}(1,1);
    upperBorder = peakBorder{iPeak,1}(2,1);
    if lowerBorder == 0 || upperBorder == 0 %case for empty peak
        peakXIC{iPeak,1} = [];
    else
        peakXIC{iPeak,1} = xic(lowerBorder:upperBorder,:);
    end
end
%reshape to original form
peakXIC = reshape(peakXIC,originalnRows,originalnCols);
end