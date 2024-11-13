function Scans = ConvertScans2MolecularMass(Scans,polarities)
% use scan polarities stored in polarities to convert MS scans from
% pseudomolecular mass to molecular mass

arguments
            Scans       (:,1) cell
            polarities  (:,1) string
end

test = height(Scans) ~= height(polarities);
switch test
    case true
        error("Size of Scans and polarities must match")
end

modifier = ones(size(Scans))*1.007825;
idx = polarities == "+";
modifier(idx) = modifier(idx)*-1;

parfor n = 1:height(Scans)
    if isempty(Scans{n,1})
        continue
    end
    Scans{n,1} = [Scans{n,1}(:,1) + modifier(n),Scans{n,1}(:,2)];
end
