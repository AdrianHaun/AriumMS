function FeatureStruct = storeInNewFeat(FeatureStruct,currentFileID,peakArray)
%% storeInNewFeat finds first empty entry in preallocated feature struct and stores peak there

firstEmpty = find(isnan(vertcat(FeatureStruct(:).retentionTime)));
firstEmpty = firstEmpty(1);

FeatureStruct(firstEmpty).peakLocations(currentFileID) = peakArray(1,1);
FeatureStruct(firstEmpty).retentionTimes(currentFileID) = peakArray(1,2);
FeatureStruct(firstEmpty).retentionTime = peakArray(1,2);
FeatureStruct(firstEmpty).peakBorders(:,currentFileID) = [peakArray(1,3);peakArray(1,4)];
FeatureStruct(firstEmpty).peakHeights(currentFileID) = peakArray(1,5);
FeatureStruct(firstEmpty).peakAreas(currentFileID) = peakArray(1,6);
FeatureStruct(firstEmpty).entropy(currentFileID) = peakArray(1,7);
FeatureStruct(firstEmpty).signal2Noise(currentFileID) = peakArray(1,8);
FeatureStruct(firstEmpty).asymmetry = peakArray(1,10);