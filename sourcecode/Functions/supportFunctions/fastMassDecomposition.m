function decompositions = fastMassDecomposition(M, epsilon)
% integer+Range+heuristik+sort

scalingFactor = 1e5;
M_scaled = int32(round(M * scalingFactor));
epsilon_scaled = int32(round(epsilon * scalingFactor));

massVector = {12.000000, 1.007825, 78.918338, 34.968853, 18.998403, 126.904468, 14.003074, 15.994915, 30.973762, 31.972071}';
massInt = cellfun(@(x) int32(round(x * scalingFactor)), massVector, 'UniformOutput', false);

elements = struct( ...
    'symbol', {"C","H","Br","Cl","F","I","N","O","P","S"}', ...
    'mass', massVector, ...
    'massInt', massInt);

nElements = numel(elements);

maxCounts = [ ...
    90, ...% C
    160,...% H
    4, ...% Br
    6, ...% Cl
    5, ...% F
    3, ...% I
    25,...% N
    30,...% O
    6, ...% P
    6  ...% S
];

[~, perm] = sort([elements.mass], 'ascend');
invPerm = zeros(size(perm));
invPerm(perm) = 1:length(perm);

elements = elements(perm);
maxCounts = maxCounts(perm);

currentCounts = zeros(1, nElements, 'int32');

% Preallocate a big matrix to hold solutions
maxSolutions = int32(1e7); % adjust if needed
solutionsMatrix = zeros(maxSolutions, nElements, 'int32');
solutionIdx = int32(0);

% Precompute max mass cumulative sums
maxMassCum = zeros(1, nElements, 'int32');
cumsumMax = int32(0);
for i = 1:nElements
    cumsumMax = cumsumMax + int32(maxCounts(i)) * elements(i).massInt;
    maxMassCum(i) = cumsumMax;
end

% Start recursion
[solutionsMatrix, solutionIdx] = iterativeDecompose(M_scaled, epsilon_scaled, elements, maxCounts, currentCounts, nElements, solutionsMatrix, solutionIdx, maxMassCum);

% Extract only the filled rows
validSolutions = solutionsMatrix(1:solutionIdx, :);

% Reorder to original element order
decompositions = validSolutions(:, invPerm);
end

function [solutionsMatrix, solutionIdx] = iterativeDecompose(targetMass, epsilon, elements, maxCounts, counts, k, solutionsMatrix, solutionIdx,maxMassCum) %#codegen

if k == 0
    if abs(double(targetMass)) <= double(epsilon)
        solutionIdx = solutionIdx + 1;
        solutionsMatrix(solutionIdx, :) = counts;
    end
    return
end

mass_k = elements(k).massInt;

max_k = min((targetMass + epsilon) / mass_k, maxCounts(k));

if max_k < 0
    return
end

minMass = int32(0);
maxMass = int32(0);
if k > 1
    maxMass = maxMassCum(k-1);
end

for count_k = 0:max_k
    residualMass = targetMass - int32(count_k)*mass_k;

    if residualMass + epsilon < minMass
        break
    end
    if residualMass - epsilon > maxMass
        continue
    end

    newCounts = counts;
    newCounts(k) = int32(count_k);

    [solutionsMatrix, solutionIdx] = iterativeDecompose(residualMass, epsilon, elements, maxCounts, newCounts, k-1, solutionsMatrix, solutionIdx, maxMassCum);
end

end
