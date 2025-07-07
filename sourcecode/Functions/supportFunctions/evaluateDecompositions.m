function evaluation = evaluateDecompositions(decompositions, queryMass)
%EVALUATEDECOMPOSITIONS Summary of this function goes here
%   Detailed explanation goes here

isValid = false(height(decompositions),1);
RDBE = zeros(height(decompositions),1);

% check formulas
parfor iDecomposition = 1:height(decompositions)
    [isValid(iDecomposition,1), RDBE(iDecomposition,1)] = validateFormula(decompositions(iDecomposition,:));
end

%remove non Valid formulas
decompositions = decompositions(isValid,:);
RDBE = RDBE(isValid);

%% build Formula string and calculate masses and deviation from searched mass

formulaStr = strings(height(decompositions),1);
formulaMass = zeros(height(decompositions),1);

parfor iDecomposition = 1:height(decompositions)
    [formulaStr(iDecomposition,1), formulaMass(iDecomposition,1)] = decompositionToFormula(decompositions(iDecomposition,:));
end

massDeviationPPM = (abs(formulaMass-queryMass)./queryMass)*10^6;

score = calculateScore(massDeviationPPM,RDBE,decompositions);

varNames = ["Formula","mass deviation [ppm]","Formula Mass","Score","double-bond equivalents"];

evaluation = table(formulaStr,massDeviationPPM,formulaMass,score,RDBE ,'VariableNames',varNames);
evaluation(evaluation.Score < 0.15,:) = [];
evaluation = sortrows(evaluation,"Score","descend");

end

function score = calculateScore(massDeviation,valuesRDBE,decompositions)

rareElementSum = zeros(height(decompositions),1);

w1 = 100;
w2 = 2;
w3 = 5;

parfor iDecmoposition = 1:height(decompositions)
    % Elemente
    counts = double(decompositions(iDecmoposition,:));

    Br = counts(3) * 3;
    Cl = counts(4) * 2;
    F = counts(5) * 3;
    I = counts(6) * 4;
    N = counts(7) * 0.5;
    P = counts(9) * 1;
    S = counts(10) * 1;

    rareElementSum(iDecmoposition) = sum(Br + Cl + F + I + P + S + N);
end
score = w1*(1-massDeviation) + w2*valuesRDBE + w3*rareElementSum;

score = 1./score;
score = score./max(score);

end