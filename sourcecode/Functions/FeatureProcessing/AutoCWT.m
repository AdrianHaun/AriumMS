function peaks = AutoCWT(DiffEIC,smoothedEIC,FilterBank)
% performs CWT and initial peak picking on ROI

%% wavelet transform
[CWT,~,~,~] = wt(FilterBank,-DiffEIC);
CWT=rescale(real(CWT),0,1);
%find initial rt
RTID=any(imextendedmax(CWT,0.2),1);
smoothedEIC(~RTID)=0;
[~,locs,~,~] = findpeaks(smoothedEIC,'WidthReference','halfheight');
clearvars smoothedEIC DiffEIC
%find initial border locations
[~,borders,~,~] = findpeaks(sum(imextendedmin(CWT,0.1)));
clearvars CWT
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