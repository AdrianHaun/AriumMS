function peakXIC = ExtractXIC(xic,peakBorder)
%reshape input to column vector
[originalnRows,originalnCols] = size(peakBorder);
peakBorder = reshape(peakBorder,[],1);
%preallocation
peakXIC = cell(size(peakBorder));
for numPeak = 1:size(peakXIC,1)
    lowerBorder = peakBorder{numPeak,1}(1,1);
    upperBorder = peakBorder{numPeak,1}(2,1);
    if lowerBorder == 0 || upperBorder == 0 %case for empty peak
        peakXIC{numPeak,1} = [];
    else
        peakXIC{numPeak,1} = xic(lowerBorder:upperBorder,:);
    end
end
%reshape to original form
peakXIC = reshape(peakXIC,originalnRows,originalnCols);
end