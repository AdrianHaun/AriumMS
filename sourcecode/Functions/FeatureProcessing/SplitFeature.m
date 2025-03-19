function featuresOut = SplitFeature(featureIn,idToSplit)

arguments
    featureIn struct
    idToSplit (1,:) double
end

if isempty(idToSplit)
    featuresOut = featureIn;
    return
end

nFiles = numel(featureIn.peakHeights);

%empty feature
emptyFeature = featureIn([],1);
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

featuresOut = repmat(emptyFeature,numel(idToSplit),1);

for n = 1:numel(idToSplit)
    featuresOut(n).mass_measured = featureIn.mass_measured;
    featuresOut(n).peakHeights(idToSplit(n)) =  featureIn.peakHeights(idToSplit(n));
    featuresOut(n).peakAreas(idToSplit(n)) =  featureIn.peakAreas(idToSplit(n));
    featuresOut(n).peakLocations(idToSplit(n)) =  featureIn.peakLocations(idToSplit(n));
    featuresOut(n).peakBorders(:,idToSplit(n)) =  featureIn.peakBorders(:,idToSplit(n));
    featuresOut(n).retentionTimes(idToSplit(n)) =  featureIn.retentionTimes(idToSplit(n));
    featuresOut(n).retentionTime(idToSplit(n)) =  featureIn.retentionTimes(idToSplit(n));
    featuresOut(n).signal2Noise(idToSplit(n)) =  featureIn.signal2Noise(idToSplit(n));
    featuresOut(n).entropy(idToSplit(n)) =  featureIn.entropy(idToSplit(n));
    featuresOut(n).XIC(idToSplit(n)) =  featureIn.XIC(idToSplit(n));
    featuresOut(n).spectrumMS1(idToSplit(n)) =  featureIn.spectrumMS1(idToSplit(n));
    featuresOut(n).spectrumMS2(idToSplit(n)) =  featureIn.spectrumMS2(idToSplit(n));
end

%remove Peaks from Input Feature
featureIn.peakHeights(idToSplit) = NaN;
featureIn.peakAreas(idToSplit) = NaN;
featureIn.peakLocations(idToSplit) = NaN;
featureIn.peakBorders(:,idToSplit) = NaN(2,1);
featureIn.retentionTimes(idToSplit) = NaN;
featureIn.signal2Noise(idToSplit) = NaN;
featureIn.entropy(idToSplit) = NaN;
featureIn.XIC(idToSplit) = cell(1);
featureIn.spectrumMS1(idToSplit) = cell(1);
featureIn.spectrumMS2(idToSplit) = cell(1);
featureIn.retentionTime = mean(featureIn.retentionTimes,"all","omitmissing");

%append
featuresOut = [featureIn;featuresOut];
