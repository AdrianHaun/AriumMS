function filteredScan = InScanIsotopeFilter(ScanCells)

filteredScan = cell(size(ScanCells));
neutronMass = 1.00866491606;
massTolerance = 0.01;

parfor n = 1:height(ScanCells)
    Scan = ScanCells{n,1};
    Scan = flipud(Scan);
    for m = 1:height(Scan)
        currentIsotope = Scan(m,1);
        id = abs(Scan(:,1)-(currentIsotope-neutronMass)) <= massTolerance;
        id = Scan(id,2) > Scan(m,2);
        if any(id)
            Scan(m,2) = NaN;
        end
    end
    Scan(isnan(Scan(:,2)),:) = [];
    filteredScan{n,1} = flipud(Scan);
end