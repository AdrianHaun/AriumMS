function EIC = buildEIC(peakCells,mz,tol,tolUnit)
%PLOTEIC Summary of this function goes here
%   Detailed explanation goes here

EIC = zeros(size(peakCells));

for n = 1:size(peakCells,1)
    scan = peakCells{n,1};
    switch tolUnit
        case "Da"
            id = abs(scan(:,1)-mz)<=tol;
        case "ppm"
            id = abs(((scan(:,1)-mz)./mz)*10^6)<=tol;
    end
    EIC(n,1) = sum(scan(id,2));
end
EIC = smoothdata(EIC,1,"gaussian",20,"omitmissing");
