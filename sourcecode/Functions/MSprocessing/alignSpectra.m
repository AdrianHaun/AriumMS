function alignedSpectra = alignSpectra(spectraCells,mode,resolution,normalize,threshold)
%% alignSpectra takes a cell array of mass spectra and aligns them along the mass axis
% 
% inputs:   spectraCells: cell array containing two column matrices, 
%                         column1: mass; column 2: intensity
% optional inputs: 
% mode: Changes the output. Default is "normal"
%       "normal": aligns spectra, output is a 1+n column matrix
%       "average": calculates mean spectra, output is a 2 column matrix.
% resolution: changes which mass delta is considered a different mass. 
%             Default is "high"
%             "high": 0.01 Da, "medium": 0.05 Da, "low": 0.1 Da
% normalize: normalizes each spectra intensities between 0 and 1, and
%            removes intensities < 1%. Default is "true
%            "true", normalize intensities
%            "false", keep original intensities
%
% output:   matrix, first column: masses, other columns: spectra intensities
%           either two column when mode = "average" or 
%           1+n when mode = "normal" with n = number of input spectra

arguments
    spectraCells    (:,1) cell
    mode            (1,1) string {mustBeMember(mode,["normal","average"])} = "normal"
    resolution      (1,1) string {mustBeMember(resolution,["low","medium","high"])} = "high"
    normalize       (1,1) string {mustBeMember(normalize,["false","true"])} = "true"
    threshold      (1,1) double {mustBeInRange(threshold,0,1)} = 0.01
end

switch resolution
    case "high" 
        mzerror = 0.01;
    case "medium"
        mzerror = 0.05;
    case "low"
        mzerror = 0.10;
end

%remove empty cells
spectraCells(cellfun(@isempty,spectraCells)) = [];

if isscalar(spectraCells) %only one spectra-> just unpack
    alignedSpectra = spectraCells{1,1};
elseif isempty(spectraCells)
    alignedSpectra = [];
else %align scans using ROI
    %build synthetic time vector
    times = 1:numel(spectraCells);
    % align
    [mzroi,MSroi,~] = ROIpeaks3(spectraCells,0,mzerror,"Da",1,times);
    if strcmp("average",mode)
        MSroi = sum(MSroi);
    end

    alignedSpectra = [mzroi;MSroi]';
end

%normalize
if strcmp(normalize,"true")
    alignedSpectra(:,2:end) = alignedSpectra(:,2:end)./max(alignedSpectra(:,2:end));
    %remove masses with rel.intensity < threshold
    id = alignedSpectra < threshold;
    alignedSpectra(id) = 0;
    id = all(alignedSpectra(:,2:end)==0,2);
    alignedSpectra(id,:) = [];
end

