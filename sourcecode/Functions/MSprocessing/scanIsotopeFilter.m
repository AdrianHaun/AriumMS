function deisotopedScan = scanIsotopeFilter(originalScan,massTolerance,tolUnit)
%% scanIsotopeFilter removes isotopes from MS scans
%
% Identifies isotopes for up to five fold charged cases.
%
% inputs: originalScan: cell array containing two column matrices, 
%                    column1: mass; column 2: intensity 
%         massTolerance: changes which mass delta is accepted to be an isotope. 
%         tolUnit: which tolerance metric to use. either "Da" (absolute) or 
%         "ppm" (relative)
%
% output: deisotopedScans: deisotoped scans in the same format as the input

arguments
    originalScan (:,1) cell
    massTolerance (1,1) double = 0.01
    tolUnit (1,1) string {mustBeMember(tolUnit,["Da","ppm"])} = "Da"
end

deisotopedScan = cell(size(originalScan));
NEUTRON_MASS = 1.00866491606;

parfor iScan = 1:height(originalScan)
    currentScan = originalScan{iScan,1};
    currentScan = flipud(currentScan);
    id = [];
    for jMass = 1:height(currentScan)
        currentIsotope = currentScan(jMass,1);
        switch tolUnit
            case "Da"
                id = abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS)) <= massTolerance |...
                    abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS/2)) <= massTolerance |...
                    abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS/3)) <= massTolerance |...
                    abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS/4)) <= massTolerance |...
                    abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS/5)) <= massTolerance;
            otherwise
                id = abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS))/(currentIsotope-NEUTRON_MASS)*1e6 <= massTolerance |...
                    abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS/2))/(currentIsotope-NEUTRON_MASS/2)*1e6 <= massTolerance |...
                    abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS/3))/(currentIsotope-NEUTRON_MASS/3)*1e6 <= massTolerance |...
                    abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS/4))/(currentIsotope-NEUTRON_MASS/4)*1e6 <= massTolerance |...
                    abs(currentScan(:,1)-(currentIsotope-NEUTRON_MASS/5))/(currentIsotope-NEUTRON_MASS/5)*1e6 <= massTolerance;
        end
        id = currentScan(id,2) > currentScan(jMass,2);
        if any(id)
            currentScan(jMass,2) = NaN;
        end
    end
    currentScan(isnan(currentScan(:,2)),:) = [];
    deisotopedScan{iScan,1} = flipud(currentScan);
end