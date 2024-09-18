function filteredScan = InScanIsotopeFilter(ScanCells,massTolerance,tolUnit)

filteredScan = cell(size(ScanCells));
neutronMass = 1.00866491606;

parfor n = 1:height(ScanCells)
    Scan = ScanCells{n,1};
    Scan = flipud(Scan);
    id = [];
    for m = 1:height(Scan)
        currentIsotope = Scan(m,1);
        switch tolUnit
            case "Da"
                id = abs(Scan(:,1)-(currentIsotope-neutronMass)) <= massTolerance;
            case "ppm"
                id = abs(Scan(:,1)-(currentIsotope-neutronMass))/(currentIsotope-neutronMass)*1e6 <= massTolerance;
        end
        id = Scan(id,2) > Scan(m,2);
        if any(id)
            Scan(m,2) = NaN;
        end
    end
    Scan(isnan(Scan(:,2)),:) = [];
    filteredScan{n,1} = flipud(Scan);
end