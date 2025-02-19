function IntResults = FinalizeIntegrationOutput(IntResults,EIC,times)
% Performs Integration of found Peaks and gathers retention times
%check for empty struct
if isempty(IntResults)
    return
end
areas=zeros(size(IntResults.peakLocation));
EIC = full(EIC);
for n=1:numel(areas)
    areas(n,1)=trapz(EIC(IntResults.peakStartLocation(n,1):IntResults.peakEndLocation(n,1)));
end
IntResults.peakArea = areas;
IntResults.peakRetentionTime = times(IntResults.peakLocation);
end