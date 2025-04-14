function eic = buildEIC(peakCells,mass,tolerance,toleranceUnit)
%% buildEIC builds extracted ion chromatogram from mass scans
%   Searches in each mass scan for an mz value within the supplied
%   tolerance and build the extracted ion chromatogram intensities, then
%   apply Gaussian smoothing.
%
% inputs: peakCells: cell array containing two column matrices, 
%                    column1: mass; column 2: intensity 
%
%         mz: mass value to search
%
%         tolerance: mass tolerance to accept; either in Da or ppm
%
%         toleranceUnit: which tolerance metric to use.
%         either "Da" (absolute) or "ppm" (relative)
%
% output: double vector containing the summed intensities of all masses
%         within tolerance of each scan

arguments
    peakCells       (:,1) cell
    mass            (1,1) double {mustBeFinite,mustBePositive}
    tolerance       (1,1) double {mustBeFinite,mustBePositive}
    toleranceUnit   (1,1) string {mustBeMember(toleranceUnit,["Da","ppm"])}
end

eic = zeros(size(peakCells));

parfor n = 1:size(peakCells,1)
    scan = peakCells{n,1};
    switch toleranceUnit
        case "Da"
            id = abs(scan(:,1)-mass) <= tolerance;
        case "ppm"
            id = abs(((scan(:,1)-mass)./mass)*10^6) <= tolerance;
    end
    eic(n,1) = sum(scan(id,2));
end

eic = smoothdata(eic,1,"gaussian",20,"omitmissing");
