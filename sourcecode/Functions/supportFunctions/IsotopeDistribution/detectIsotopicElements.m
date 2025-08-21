function elementHits = detectIsotopicElements(isoMasses, isoIntensities, chargeState)
%DETECTISOTOPICELEMENTS Detect Cl, Br, S elements including overlapping contributions.
%
% Inputs:
%   isoMasses      - vector of isotope masses (Da)
%   isoIntensities - vector of normalized intensities (0-1)
%   chargeState    - integer charge
%
% Output:
%   elementHits - struct with estimated counts and confidence

    % If empty input
    if isempty(isoMasses) || isempty(isoIntensities)
        elementHits = struct('Cl', struct('count', 0, 'confident', false), ...
                             'Br', struct('count', 0, 'confident', false), ...
                             'S',  struct('count', 0, 'confident', false));
        return;
    end

    % Expected spacing per charge
    spacing = 1 / chargeState;

    % Find +2 Da peak
    baseMass = isoMasses(1);
    targetMass = baseMass + 2 * spacing;
    deltaMasses = abs(isoMasses - targetMass);
    [minDelta, idx2Da] = min(deltaMasses);

    if minDelta > 0.02
        I2Da = 0;
    else
        I2Da = isoIntensities(idx2Da);
    end

    I0 = isoIntensities(1);

    ratio = I2Da / I0;

    % Basis contributions per atom
    rCl = 0.33;
    rBr = 1.0;
    rS = 0.042;

    % Solve: [rCl rBr rS] * [nCl; nBr; nS] = ratio
    A = [rCl rBr rS];

    % Least squares estimate with non-negative constraint
    nAtoms = lsqnonneg(A, ratio);

    % Round counts to nearest integer
    estCl = round(nAtoms(1));
    estBr = round(nAtoms(2));
    estS  = round(nAtoms(3));

    % Recompute expected ratio
    expectedRatio = estCl*rCl + estBr*rBr + estS*rS;

    % Confidence if difference < threshold
    confThreshold = 0.5;
    conf = abs(expectedRatio - ratio) < confThreshold;

    % Output struct
    elementHits = struct();

    elementHits.Cl.count = estCl;
    elementHits.Cl.confident = conf && estCl>0;

    elementHits.Br.count = estBr;
    elementHits.Br.confident = conf && estBr>0;

    elementHits.S.count = estS;
    elementHits.S.confident = conf && estS>0;
end
