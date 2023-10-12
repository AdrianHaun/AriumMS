function peaks = AutoCWT(DiffEIC,smoothedEIC,ScanFrequency,minWidth,maxWidth)
% performs CWT and initial peak picking on ROI
FilterBank = cwtfilterbank("SignalLength",numel(DiffEIC),"WaveletParameters",[3 4],"VoicesPerOctave",8,"SamplingPeriod",seconds(ScanFrequency),"PeriodLimits",[seconds(minWidth) seconds(maxWidth)]);% prepare wavelet filterbank
%% wavelet transform
[CWT,~,~,~] = wt(FilterBank,-DiffEIC);
CWT=rescale(real(CWT),0,1);
%find initial rt
RTID=any(imextendedmax(CWT,0.2),1);
isolatedPeak=smoothedEIC;
isolatedPeak(~RTID)=0;
[~,locs,~,~] = findpeaks(isolatedPeak,'WidthReference','halfheight');
%find initial border locations
[~,borders,~,~] = findpeaks(sum(imextendedmin(CWT,0.1)));
peaks=zeros(numel(locs),4);
peaks(:,1)=locs';
%sort borders to RT
for n=1:size(peaks,1)
    [~,idx]=mink(borders-locs(n),2,'ComparisonMethod','abs');
    if ~isempty(borders(idx))==true
        peaks(n,2:3)=sort(borders(idx),'ascend');
    end
end
idx=peaks(:,2)==0;
peaks(idx,:)=[];
end