function FeatureStruct = storeInNewFeat(FeatureStruct,currentFileID,peakArray)
%% storeInNewFeat finds first empty entry in preallocated feature struct and stores peak there

firstEmpty = find(isnan(vertcat(FeatureStruct(:).rt)));
firstEmpty = firstEmpty(1);

FeatureStruct(firstEmpty).peakLocation(currentFileID) = peakArray(1,1);
FeatureStruct(firstEmpty).peakRTs(currentFileID) = peakArray(1,2);
FeatureStruct(firstEmpty).rt = peakArray(1,2);
FeatureStruct(firstEmpty).peakBorder(:,currentFileID) = [peakArray(1,3);peakArray(1,4)];
FeatureStruct(firstEmpty).Int(currentFileID) = peakArray(1,5);
FeatureStruct(firstEmpty).Area(currentFileID) = peakArray(1,6);
FeatureStruct(firstEmpty).entropy(currentFileID) = peakArray(1,7);
FeatureStruct(firstEmpty).signal2Noise(currentFileID) = peakArray(1,8);
FeatureStruct(firstEmpty).asymmetry = peakArray(1,10);