function [formulaStr, totalMass] = decompositionToFormula(counts)
% counts: Vektor der Elementanzahlen
% Rückgabe:
%   formulaStr : z.B. "C6H12O6"
%   totalMass  : berechnete Masse der Summenformel

counts = double(counts);

% Fest vorgegebene Elemente
elements = struct( ...
    'symbol', {"C","H","Br","Cl","F","I","N","O","P","S"}', ...
    'mass', { ...
        12.000000, ...
        1.007825, ...
        78.918338, ...
        34.968853, ...
        18.998403, ...
        126.904468, ...
        14.003074, ...
        15.994915, ...
        30.973762, ...
        31.972071 ...
    }');

formulaStr = '';
totalMass = 0;

for i = 1:numel(elements)
    n = counts(i);
    if n > 0
        % Summenformel-Text
        if n == 1
            formulaStr = [formulaStr, elements(i).symbol];
        else
            formulaStr = [formulaStr, elements(i).symbol, num2str(n)];
        end
        % Masse aufsummieren
        totalMass = totalMass + n * elements(i).mass;
    end
end

%process formula string
formulaStr = join(formulaStr);
formulaStr = erase(formulaStr," ");
end
