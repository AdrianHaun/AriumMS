function [cleanedDataMS1,varargout] = cleanRawProfileScans(rawDataMS1,varargin)
%% cleanRawProfileScans removes noise and centroides scans
% DESCRIPTIVE TEXT


cleanedDataMS1 = struct('profileDataMS1',[],...
    'timeDataMS1',[],...
    'polarityMS1',[]);

%clean MS1
%remove empty scans
emptyScans = cellfun(@isempty, rawDataMS1.profileDataMS1);
cleanedDataMS1.profileDataMS1(emptyScans,:) = [];
cleanedDataMS1.timeDataMS1(emptyScans,:) = [];
cleanedDataMS1.polarityMS1(emptyScans,:) = [];
% set noise to 0
cleanedDataMS1.profileDataMS1 = denoiseScans(cleanedDataMS1.profileDataMS1);

%% clean MS2
if nargin == 2
    cleanedDataMS2 = struct('profileDataMS2',[],...
        'timeDataMS2',[],...
        'polarityMS2',[],...
        'precursorMass',[],...
        'fragmentationEnergy',[],...
        'fragmentationType',[]);

    rawDataMS2 = varargin{1};
    % remove noise
    denoisedScans = denoiseScans(noisyScans)
    % remove possible empty scans in ms2
    emptyScans = cellfun(@isempty, rawDataMS2.profileDataMS2);
    cleanedDataMS2.profileDataMS2(emptyScans,:) = [];
    cleanedDataMS2.timeDataMS2(emptyScans,:) = [];
    cleanedDataMS2.polarityMS2(emptyScans,:) = [];
    cleanedDataMS2.precursorMass(emptyScans,:) = [];
    cleanedDataMS2.fragmentationEnergy(emptyScans,:) = [];
    cleanedDataMS2.fragmentationType(emptyScans,:) = [];
    % normalize

    varargout{1} = cleanedDataMS2;
end
