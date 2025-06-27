function formulas = fastMassDecomposition(targetMass, tolerance)
% FASTMASSDECOMPOSITION Effiziente Summenformel-Suche nach Rojas-Chertó et al.
%
%   targetMass - gewünschte monoisotopische Masse (double)
%   tolerance  - Toleranz in Dalton (double)
%   elements   - struct array mit Feldern 'symbol' und 'mass'
%   maxCounts  - maximale Anzahl pro Element (array)

% %element search space
% elements = struct( ...
%     'symbol', {'H','C','N','O','P','S','Cl','Br'}, ...
%     'mass', [1.007825;12;14.003074;15.994915;30.973762;31.972071;34.968853;78.918227;126.904473]);
% % element count search space
% maxCounts = [40,20,5,15,3,2,1,1];

%element search space
elements = struct( ...
    'symbol', {"C","H","N","O"}', ...
    'mass', {12,1.007825,14.003074,15.994915}');
% element count search space
maxCounts = [10,10,2,5];

scaleFactor = 1e6;
nElements = length(elements);

% Diskrete Zielmasse
massInt = round(targetMass * scaleFactor);
tolInt  = round(tolerance  * scaleFactor);
targetMin = massInt - tolInt;
targetMax = massInt + tolInt;

% Diskrete Elementmassen
elementMassesInt = round([elements.mass] * scaleFactor);

% Elemente optional nach Masse sortieren (schnelleres Pruning)
[~, sortIdx] = sort(elementMassesInt);
elementMassesInt = elementMassesInt(sortIdx);
elements = elements(sortIdx);
maxCounts = maxCounts(sortIdx);

% Initialisiere DP: jede Zelle ist eine Matrix der Kombinationen

% Speicherbereich aufteilen
dpBelow = cell(targetMin, 1);
dp = cell(targetMax - targetMin + 1, 1);

% Startkombination
startCombs = zeros(1, nElements);

% DP Iteration (äußere Schleife bleibt seriell)
for m = 0:targetMax
    if m==0
        currCombs = startCombs;
    elseif m<targetMin
        if isempty(dpBelow{m}), continue; end
        currCombs = dpBelow{m};
    else
        idxCurr = m - targetMin + 1;
        if isempty(dp{idxCurr}), continue; end
        currCombs = dp{idxCurr};
    end

    numCombs = size(currCombs, 1);
    newCombsAll = cell(1, numCombs);
    newMassesAll = cell(1, numCombs);

    % PARALLEL: Bearbeite jede Kombination gleichzeitig
    parfor i = 1:numCombs
        comb = currCombs(i,:);
        localNewCombs = [];
        localNewMasses = [];

        for e = 1:nElements
            if comb(e) < maxCounts(e)
                newComb = comb;
                newComb(e) = newComb(e) + 1;
                newMass = m + elementMassesInt(e);

                if newMass > targetMax
                    continue;
                end

                massLeft = sum((maxCounts - newComb) .* elementMassesInt);
                if newMass + massLeft < targetMin
                    continue;
                end

                localNewCombs = [localNewCombs; newComb];
                localNewMasses = [localNewMasses; newMass];
            end
        end
        newCombsAll{i} = localNewCombs;
        newMassesAll{i} = localNewMasses;
    end

    % Seriell: Ergebnisse einsortieren
    for i = 1:numCombs
        combs_i = newCombsAll{i};
        masses_i = newMassesAll{i};

        for j = 1:size(combs_i,1)
            newComb = combs_i(j,:);
            newMass = masses_i(j);

            if newMass < targetMin
                if isempty(dpBelow{newMass})
                    dpBelow{newMass} = newComb;
                else
                    if ~ismember(newComb, dpBelow{newMass}, 'rows')
                        dpBelow{newMass} = [dpBelow{newMass}; newComb];
                    end
                end
            else
                idxNew = newMass - targetMin + 1;
                if isempty(dp{idxNew})
                    dp{idxNew} = newComb;
                else
                    if ~ismember(newComb, dp{idxNew}, 'rows')
                        dp{idxNew} = [dp{idxNew}; newComb];
                    end
                end
            end
        end
    end
end



% Ergebnisse extrahieren
formulas = {};
for m = targetMin:targetMax
    idx = m - targetMin + 1;
    if isempty(dp{idx})
        continue
    end
    for i = 1:size(dp{idx},1)
        c = dp{idx}(i,:);
        c_unsorted = zeros(1, nElements);
        c_unsorted(sortIdx) = c;

        formula = '';
        for e = 1:nElements
            count = c_unsorted(e);
            if count>0
                formula = [formula, elements(sortIdx(e)).symbol];
                if count>1
                    formula = [formula, num2str(count)];
                end
            end
        end
        formulas{end+1} = formula;
    end
end
%unpack and join
formulas = vertcat(formulas{:});
formulas(:,1) = [];
formulas = join(formulas,"");

% Doppelte Summenformeln entfernen
formulas = unique(formulas);

end
