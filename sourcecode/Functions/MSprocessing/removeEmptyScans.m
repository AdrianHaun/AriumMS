function [cleanedDataMS1,varargout] = removeEmptyScans(rawDataMS1,varargin)
%% removeEmptyScans removes empty scans from MS-Data struct and compresses profileData by removing entries without intensity
% DESCRIPTIVE TEXT

cleanedDataMS1 = rawDataMS1;
%remove entries with 0 intensity
profileDataMS1 = cleanedDataMS1.profileDataMS1;
parfor iScan = 1:height(profileDataMS1)
    toRemove = profileDataMS1{iScan,1}(:,2) == 0;
    profileDataMS1{iScan,1}(toRemove,:) = []
end
cleanedDataMS1.profileDataMS1 = profileDataMS1;

%remove empty scans
toRemove = cellfun(@isempty, cleanedDataMS1.profileDataMS1) | cellfun(@isempty, cleanedDataMS1.centroidDataMS1);
cleanedDataMS1.profileDataMS1(toRemove,:) = [];
cleanedDataMS1.centroidDataMS1(toRemove,:) = [];
cleanedDataMS1.timeDataMS1(toRemove,:) = [];
if ~isempty(cleanedDataMS1.polarityMS1)
    cleanedDataMS1.polarityMS1(toRemove,:) = [];
end


%% clean MS2
if nargin == 2
    rawDataMS2 = varargin{1};
    cleanedDataMS2 = rawDataMS2;
    emptyScans = cellfun(@isempty, cleanedDataMS2.profileDataMS2);
    cleanedDataMS2.profileDataMS2(emptyScans,:) = [];
    cleanedDataMS2.timeDataMS2(emptyScans,:) = [];
    cleanedDataMS2.polarityMS2(emptyScans,:) = [];
    cleanedDataMS2.precursorMass(emptyScans,:) = [];
    cleanedDataMS2.fragmentationEnergy(emptyScans,:) = [];
    cleanedDataMS2.fragmentationType(emptyScans,:) = [];
    varargout{1} = cleanedDataMS2;
end
