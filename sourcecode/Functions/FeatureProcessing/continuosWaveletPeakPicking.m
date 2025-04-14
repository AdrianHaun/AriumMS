function initialPeaks = continuosWaveletPeakPicking(derivativeChromatogram,smoothedChromatogram,cwtFilterBank)
%% continuosWaveletPeakPicking performs CWT for initial peak picking on Chromatograms
% Performs continuos wavelet transform (CWT) on the second derivative of an
% extracted ion chromatogramm, to find initial peak locations and peak
% borders.
%
% inputs:
% derivativeChromatogram: column vector, second derivative of a
%                           chromatogram with peaks
% smoothedChromatogram: column vector, smoothed variant of the chromatogram
%
% cwtFilterBank: CWT options defined by the cwtfilterbank function
%
% output: (n,3) double matrix, column 1 contains peak location indices,
%         column 2 lower peak borders, column 3 upper peak borders

arguments
    derivativeChromatogram  (:,1) double
    smoothedChromatogram    (:,1) double {mustBeSameHeight(derivativeChromatogram,smoothedChromatogram)}
    cwtFilterBank           (1,1) cwtfilterbank
end

% continuos wavelet transform
waveletTransformMatrix = wt(cwtFilterBank,-derivativeChromatogram);
waveletTransformMatrix = rescale(real(waveletTransformMatrix),0,1);


%find initial peak location
retentionTimeIndexArray = any(imextendedmax(waveletTransformMatrix,0.2),1);
smoothedChromatogram(~retentionTimeIndexArray) = 0;
[~,peakLocations,~,~] = findpeaks(smoothedChromatogram,'WidthReference','halfheight');

%find initial border locations
[~,peakBorders,~,~] = findpeaks(sum(imextendedmin(waveletTransformMatrix,0.1)));

clearvars waveletTransformMatrix smoothedEIC derivativeEIC

initialPeaks = zeros(numel(peakLocations),3);
initialPeaks(:,1) = peakLocations';

%sort borders to RT
for iPeak = 1:size(initialPeaks,1)
    [~,idx] = mink(peakBorders-peakLocations(iPeak),2,'ComparisonMethod','abs');
    if ~isempty(peakBorders(idx))
        initialPeaks(iPeak,2:3) = sort(peakBorders(idx),'ascend');
    end
end
%remove peaks with missing borders
initialPeaks(initialPeaks(:,2)== 0 | initialPeaks(:,3)== 0,:) = [];


%% Validation function
function mustBeSameHeight(vectorA,vectorB)
% Test for equal height
if ~isequal(height(vectorA),height(vectorB))
    eid = 'Size:notEqual';
    msg = 'Height of first input must equal size of second input.';
    error(eid,msg)
end