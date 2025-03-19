function alignedSpectra = AlignSpectra(spectraCells,mode,resolution,rescale)

arguments
    spectraCells    (:,1) cell
    mode            (1,1) string {mustBeMember(mode,["normal","average"])} = "normal"
    resolution      (1,1) string {mustBeMember(resolution,["low","medium","high"])} = "high"
    rescale            (1,1) string {mustBeMember(rescale,["false","true"])} = "true"
end

if isscalar(spectraCells) %only one spectra-> just unpack, rescale and clean
    alignedSpectra = spectraCells{1,1};

else %align
    switch resolution
        case "high"
            mzerror = 0.01;
        case "medium"
            mzerror = 0.05;
        case "low"
            mzerror = 0.10;
    end
    
    errorUnit = "Da";

    %align scans using ROI
    times = 1:numel(spectraCells);
    [mzroi,MSroi,~] = ROIpeaks3(spectraCells,0,mzerror,errorUnit,1,times);

    if strcmp("average",mode)
        MSroi = mean(MSroi);
    end

    alignedSpectra = [mzroi;MSroi]';
end

%rescale
if strcmp(rescale,"true")
    alignedSpectra(:,2:end) = alignedSpectra(:,2:end)./max( alignedSpectra(:,2:end));
end

%remove masses with rel.intensity < 1%
id = alignedSpectra < 0.01;
alignedSpectra(id) = 0;
id = all(alignedSpectra(:,2:end)==0,2);
alignedSpectra(id,:) = [];
end


