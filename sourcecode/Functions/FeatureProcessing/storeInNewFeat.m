function currentFeatureStruct = storeInNewFeat(currentFeatureStruct,currentFileID,peakData)
% finds first empty entry in preallocated feature struct and stores peak
% there

firstEmpty = find(isnan(vertcat(currentFeatureStruct(:).retentionTime)));
firstEmpty = firstEmpty(1);
currentFeatureStruct(firstEmpty).peakLocations(currentFileID) = peakData(1,1);
currentFeatureStruct(firstEmpty).retentionTimes(currentFileID) = peakData(1,2);
currentFeatureStruct(firstEmpty).retentionTime = peakData(1,2);
currentFeatureStruct(firstEmpty).peakBorders(:,currentFileID) = [peakData(1,3);peakData(1,4)];
currentFeatureStruct(firstEmpty).peakHeights(currentFileID) = peakData(1,5);
currentFeatureStruct(firstEmpty).peakAreas(currentFileID) = peakData(1,6);
currentFeatureStruct(firstEmpty).entropy(currentFileID) = peakData(1,7);
currentFeatureStruct(firstEmpty).signal2Noise(currentFileID) = peakData(1,8);
currentFeatureStruct(firstEmpty).asymmetry = peakData(1,10);