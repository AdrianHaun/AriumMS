function compositScore = calculateCompositeScore(spectrumA,spectrumB)
%% calculateCompositScore calculates the composite score between two spectra
%
% inputs: spectrumA, spectrumB: spectra to compare. Specified as
%                               double column vector of intensity values
%                               with the same mass axis
% output: composit score similarity measure. (1,1) double, between 0 and
% 999
%
arguments
    spectrumA (:,1) {mustBeNumeric,mustBeReal}
    spectrumB (:,1) {mustBeNumeric,mustBeReal,mustBeEqualSize(spectrumA,spectrumB)}
end

SCORE_MAXIMUM = 999;

nZeroSpectrumA = sum(spectrumA ~= 0); %number of non zero elements in spectrumA
nZeroBothSpectra = sum(spectrumA ~= 0 & spectrumB ~= 0); %number of non zero elements in both vectors

% calculate cosine similarity
cosineSimilarity = dot(spectrumA,spectrumB)^2/(sum(spectrumA.^2)*sum(spectrumB.^2));

%remove zero elements
indexToRemove = spectrumA == 0 | spectrumB == 0;

spectrumA(indexToRemove) = [];
spectrumB(indexToRemove) = [];

% shift vectors by 1
downShiftedA = [0;spectrumA];
downShiftedB = [0;spectrumB];
spectrumA = [spectrumA;0];
spectrumB = [spectrumB;0];

% calculate correction factor
composit = (spectrumA./downShiftedA.*downShiftedB./spectrumB);
composit(~isfinite(composit)) = [];
composit(composit>1) = composit(composit>1).^-1;

correctionFactor = 1/nZeroBothSpectra * sum(composit);

%% calculate composite score
compositScore = SCORE_MAXIMUM*(nZeroSpectrumA*cosineSimilarity + nZeroBothSpectra*correctionFactor)/(nZeroSpectrumA + nZeroBothSpectra);
compositScore = round(compositScore);

% Custom validation function
function mustBeEqualSize(a,b)
% Test for equal size
if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    error(eid,msg)
end