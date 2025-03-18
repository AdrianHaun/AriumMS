function alignedSpectra = AlignSpectra(spectraCells,mode)

arguments
    spectraCells (:,1) cell
    mode (1,1) string {mustBeMember(mode,["normal","average"])} = "normal"
end

if isscalar(spectraCells) %only one spectra-> just unpack
    alignedSpectra = spectraCells{1,1};

else %align

    mzerror = 0.10;
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
alignedSpectra(:,2:end) = alignedSpectra(:,2:end)./max( alignedSpectra(:,2:end));
end


