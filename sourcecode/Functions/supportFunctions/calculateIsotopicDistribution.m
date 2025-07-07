function massDistribution = calculateIsotopicDistribution(decompositions)

arguments
    decompositions (:,10)
end

decompositions = double(decompositions);
massDistribution = cell(height(decompositions),1);

parfor iDecomposition = 1:height(decompositions)
    currentDecomposition = decompositions(iDecomposition,:);
    formulaStruct = struct();
    formulaStruct.C = currentDecomposition(1);
    formulaStruct.H  = currentDecomposition(2);
    formulaStruct.Br = currentDecomposition(3);
    formulaStruct.Cl = currentDecomposition(4);
    formulaStruct.F = currentDecomposition(5);
    formulaStruct.I = currentDecomposition(6);
    formulaStruct.N = currentDecomposition(7);
    formulaStruct.O = currentDecomposition(8);
    formulaStruct.P = currentDecomposition(9);
    formulaStruct.S = currentDecomposition(10);
    massDistribution{iDecomposition,1} = isotopicdist(formulaStruct);
end