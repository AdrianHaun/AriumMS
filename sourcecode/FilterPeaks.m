function IntResults = FilterPeaks(IntResults,MinPWDataPoints,MaxPWDataPoints,maxSN,Noise,EIC)
% Filters identified peaks from AutoCWT
%check empty input
if isempty(IntResults.peakLocation)
    return
end
%preallocate indexarray
idx = false(size(IntResults.peakLocation));

%% Peak filter
%remove duplicate peaks
[~,id] = unique([IntResults.peakLocation,IntResults.peakStartLocation,IntResults.peakEndLocation,IntResults.peakHeight,],'rows');
idx = idx | id;

%remove peaks with wrong boundaries
id = IntResults.peakStartLocation>=IntResults.peakEndLocation;
idx = idx | id;

%remove peaks with height = 0
id = IntResults.peakHeight == 0;
idx = idx | id;

%remove peaks with bad Peak asymmetry

symmetry = (IntResults.peakEndLocation - IntResults.peakLocation)./(IntResults.peakLocation - IntResults.peakStartLocation);
id = symmetry<0.3 | symmetry>3;
idx = idx | id;

%less than minimum peak width
id = IntResults.peakEndLocation-IntResults.peakStartLocation < MinPWDataPoints;
IntResults.minWidthFilteres=sum(id);
idx = idx | id;

%more than maximum peak width
id=IntResults.peakEndLocation - IntResults.peakStartLocation > MaxPWDataPoints;
IntResults.maxWidthFilteres=sum(id);
idx = idx | id;

%S/N peak rejection
IntResults.SN = IntResults.PeakHeight ./ Noise;
id=SN<maxSN;
IntResults.SNFiltered = sum(id);
idx = idx | id;

% remove identified peaks
IntResults.peakLocation(idx) = [];
IntResults.peakStartLocation(idx) = [];
IntResults.peakEndLocation(idx) = [];
IntResults.peakHeight(idx) = [];
IntResults.SN(idx) = [];

end