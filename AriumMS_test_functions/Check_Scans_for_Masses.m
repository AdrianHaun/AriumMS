% find mass and plot

queryMass = 118.026615000000;

%profile
ScanCells = ScansProfile;
figure
hold on
title(queryMass + " profile Scans")
for n = 1:6
    scan = ScanCells{n,1};
    EIC = buildEIC(scan,queryMass,0.02,"Da");
    plot(ScanCells{n,2},EIC)
end

%centroid
ScanCells = ScansPreIsotope;
figure
hold on
title(queryMass + " centroid Scans")
for n = 1:6
    scan = ScanCells{n,1};
    EIC = buildEIC(scan,queryMass,0.02,"Da");
    plot(ScanCells{n,2},EIC)
end

%post isotope filter
ScanCells = ScansPostIsotope;
figure
hold on
title(queryMass + " isotope filtered scans")
for n = 1:6
    scan = ScanCells{n,1};
    EIC = buildEIC(scan,queryMass,0.02,"Da");
    plot(ScanCells{n,2},EIC)
end
