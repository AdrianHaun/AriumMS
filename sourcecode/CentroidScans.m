function CentroidedScans = CentroidScans(ScanData)

CentroidedScans = cell(size(ScanData));
parfor n = 1:height(CentroidedScans)
    Scan = ScanData{n,1};
    [~,edges] = histcounts(Scan(:,2));
    idx = Scan(:,2) <= edges(2)/2;
    Scan(idx,2) = 0;
    %find mass peaks and width
    [maxInt,maxMZ,width] = findpeaks(Scan(:,2),Scan(:,1),'WidthReference','halfheight');
    %expand maxInt
    maxInt = [zeros(height(maxInt),1),maxInt];

    %compute weighted mean of mz and intensity 
    for p = 1:numel(width)
        currentMZ = maxMZ(p);
        currentWidth = width(p);
        id = abs(Scan(:,1)-currentMZ)<=currentWidth;
        vec = Scan(id,:);
        mz = mean(vec(:,1),Weights = vec(:,2)/max(vec(:,2)));
        maxInt(p,1) = mz; 
    end
    CentroidedScans{n,1} = maxInt;
end