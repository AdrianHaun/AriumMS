function FeaturesOut = splitFeature(FeatureIn,idToSplit)
%% splitFeature stores data of peaks identified to be wrongly placed in a feature in a new one, removing it from the current feature
% inputs: 
%       FeatureIn: single feature struct with peaks to remove
%       idToSplit: index to peaks that need to be stored in a different feature
% output: 
%       FeaturesOut: split features as struct array

arguments
    FeatureIn struct
    idToSplit (1,:) double
end

if isempty(idToSplit)
    FeaturesOut = FeatureIn;
    return
end

nFiles = numel(FeatureIn.peakHeights);

%empty feature
emptyFeature = FeatureIn([],1);
emptyFeature(1).adductType = strings(1);
emptyFeature(1).formula = strings(1);
emptyFeature(1).peakHeights = NaN(1,nFiles);
emptyFeature(1).peakAreas = NaN(1,nFiles);
emptyFeature(1).peakLocations = NaN(1,nFiles);
emptyFeature(1).peakBorders = NaN(2,nFiles);
emptyFeature(1).retentionTimes = NaN(1,nFiles);
emptyFeature(1).signal2Noise = NaN(1,nFiles);
emptyFeature(1).entropy = NaN(1,nFiles);
emptyFeature(1).XIC = cell(1,nFiles);
emptyFeature(1).spectrumMS1 = cell(1,nFiles);
emptyFeature(1).spectrumMS2 = cell(1,nFiles);
emptyFeature(1).isotopePattern = cell(1,nFiles);
emptyFeature(1).chargeState = NaN(1,nFiles);

FeaturesOut = repmat(emptyFeature,numel(idToSplit),1);

for iSplitIndex = 1:numel(idToSplit)
    FeaturesOut(iSplitIndex).mass_measured = FeatureIn.mass_measured;
    FeaturesOut(iSplitIndex).peakHeights(idToSplit(iSplitIndex)) =  FeatureIn.peakHeights(idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).peakAreas(idToSplit(iSplitIndex)) =  FeatureIn.peakAreas(idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).peakLocations(idToSplit(iSplitIndex)) =  FeatureIn.peakLocations(idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).peakBorders(:,idToSplit(iSplitIndex)) =  FeatureIn.peakBorders(:,idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).retentionTimes(idToSplit(iSplitIndex)) =  FeatureIn.retentionTimes(idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).retentionTime(idToSplit(iSplitIndex)) =  FeatureIn.retentionTimes(idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).signal2Noise(idToSplit(iSplitIndex)) =  FeatureIn.signal2Noise(idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).entropy(idToSplit(iSplitIndex)) =  FeatureIn.entropy(idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).XIC(idToSplit(iSplitIndex)) =  FeatureIn.XIC(idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).spectrumMS1(idToSplit(iSplitIndex)) =  FeatureIn.spectrumMS1(idToSplit(iSplitIndex));
    if ~isempty(FeatureIn.spectrumMS2)
        FeaturesOut(iSplitIndex).spectrumMS2(idToSplit(iSplitIndex)) =  FeatureIn.spectrumMS2(idToSplit(iSplitIndex));
    end
    FeaturesOut(iSplitIndex).isotopePattern(idToSplit(iSplitIndex)) =  FeatureIn.isotopePattern(idToSplit(iSplitIndex));
    FeaturesOut(iSplitIndex).chargeState(idToSplit(iSplitIndex)) =  FeatureIn.chargeState(idToSplit(iSplitIndex));
end

%remove Peaks from Input Feature
FeatureIn.peakHeights(idToSplit) = NaN;
FeatureIn.peakAreas(idToSplit) = NaN;
FeatureIn.peakLocations(idToSplit) = NaN;
FeatureIn.peakBorders(:,idToSplit) = NaN(2,numel(idToSplit));
FeatureIn.retentionTimes(idToSplit) = NaN;
FeatureIn.signal2Noise(idToSplit) = NaN;
FeatureIn.entropy(idToSplit) = NaN;
FeatureIn.XIC(idToSplit) = cell(1);
FeatureIn.spectrumMS1(idToSplit) = cell(1);
FeatureIn.spectrumMS2(idToSplit) = cell(1);
FeatureIn.isotopePattern(idToSplit) = cell(1);
FeatureIn.chargeState(idToSplit) = 0;
FeatureIn.retentionTime = mean(FeatureIn.retentionTimes,"all","omitmissing");

%append
FeaturesOut = [FeatureIn;FeaturesOut];
