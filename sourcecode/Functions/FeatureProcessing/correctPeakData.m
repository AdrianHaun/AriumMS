function correctedPeakArray = correctPeakData(initialPeakArray,originalChromatogram,smoothedChromatogram,minPeakWidth)
%% correctPeakData updates inital peak location and borders
% Performs friction border correction and moving standard deviation Border
% correction on initial peak borders from continuosWaveletPeakPicking. Then
% updates the peak location and stores the peak height
%
% inputs: 
% initialPeakArray: output from continuosWaveletPeakPicking. List of initial 
%                   peak data. Column 1: peak location (index), column 2:       
%                   lower peak border (index), upper peak border (index)
%
% originalChromatogram: Chromatogram used to pick peaks. Column vector with 
%                   intensity values
%
% smoothedChromatogram: same chromatogram as originalChromatogram but with
%                       applied smoothing
%
% output: same layout as initialPeakArray, with column 4 added 
%         (peak height) and corrected entries

arguments
    initialPeakArray        (:,3) double
    originalChromatogram    (:,1) double
    smoothedChromatogram    (:,1) double {mustBeEqualSize(originalChromatogram,smoothedChromatogram)}
    minPeakWidth            (1,1) double
end
INTENSITY_THRESHOLD_PERCENT =  0.0001;

% initially filter duplicates and bad borders
initialPeakArray = unique(initialPeakArray,"rows","stable");
hasBadBorder = initialPeakArray(:,3)-initialPeakArray(:,2) <= minPeakWidth;
initialPeakArray(hasBadBorder,:) = [];

%preallocate output
correctedPeakArray = ones(height(initialPeakArray),4);

%scale chromatogram
smoothedChromatogram = smoothedChromatogram./max(smoothedChromatogram);

%normalize threshhold
intensityThresholdAbsolute = (max(smoothedChromatogram)-min(smoothedChromatogram))*INTENSITY_THRESHOLD_PERCENT;

%% friction border correction
for iPeak = 1:size(initialPeakArray,1)
    upperBorder = initialPeakArray(iPeak,3);
    while upperBorder < numel(smoothedChromatogram)-1 %upper bond
        if smoothedChromatogram(upperBorder)-smoothedChromatogram(upperBorder+1) > intensityThresholdAbsolute
            upperBorder = upperBorder+1;
        else
            initialPeakArray(iPeak,3) = upperBorder;
            break
        end
    end
    lowerBorder = initialPeakArray(iPeak,2);
    while lowerBorder > 1 %lower bond
        if smoothedChromatogram(lowerBorder)-smoothedChromatogram(lowerBorder-1) > intensityThresholdAbsolute
            lowerBorder = lowerBorder-1;
        else
            initialPeakArray(iPeak,2) = lowerBorder;
            break
        end
    end
end

noise = originalChromatogram-smoothedChromatogram;
clearvars smoothedChromatogram

%moving STD correction to remove errors from smoothed peaks
for iPeak = 1:size(initialPeakArray,1)

    %upper bond
    upperBorder = initialPeakArray(iPeak,3);
    while upperBorder > 1 && upperBorder < numel(originalChromatogram)-8 
        if abs(originalChromatogram(upperBorder-1)-mean(originalChromatogram(upperBorder-1:upperBorder+7,1))) <= std(noise(upperBorder-1:upperBorder+7,1)) || round(originalChromatogram(upperBorder-1,1)-originalChromatogram(upperBorder,1)) == 0
            upperBorder = upperBorder-1;
        else
            correctedPeakArray(iPeak,3)=upperBorder;
            break
        end
    end

    %lower brder
    lowerBorder = initialPeakArray(iPeak,2);
    while lowerBorder <= numel(originalChromatogram)-1 && lowerBorder > 8 
        if abs(originalChromatogram(lowerBorder+1,1)-mean(originalChromatogram(lowerBorder-7:lowerBorder+1,1))) <= std(noise(lowerBorder-7:lowerBorder+1,1)) || round(originalChromatogram(lowerBorder,1)-originalChromatogram(lowerBorder+1,1)) == 0
            lowerBorder = lowerBorder+1;
        else
            correctedPeakArray(iPeak,2) = lowerBorder;
            break
        end
    end
end
%get final peak location and height
for iPeak = 1:size(initialPeakArray,1)
    temporaryChromatogam = originalChromatogram;
    temporaryChromatogam(1:correctedPeakArray(iPeak,2)-1) = 0;
    temporaryChromatogam(correctedPeakArray(iPeak,3)+1:end) = 0;
    [correctedPeakArray(iPeak,4),correctedPeakArray(iPeak,1)] = max(temporaryChromatogam);
end

% Custom validation function
function mustBeEqualSize(a,b)
% Test for equal size
if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    error(eid,msg)
end