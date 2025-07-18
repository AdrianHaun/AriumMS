function [isValid, RDBE] = validateFormula(counts)
% uses heuristics described by Tobias Kind & Oliver Fiehn
% Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry
% BMC Bioinformatics. 2007 Mar 27;8:105. doi: 10.1186/1471-2105-8-105
% to 
% counts: vector of atom counts in this order: [C,H,Br,Cl,F,I,N,O,P,S]
%
% output:
%   isValid (true/false)
%   RDBE (calculated ring double bond equivalents)

if numel(counts) ~= 10
    error('counts must have 10 elements: C,H,Br,Cl,F,I,N,O,P,S');
end
counts = double(counts);

% Elemente
C = counts(1);
H = counts(2);
Br = counts(3);
Cl = counts(4);
F = counts(5);
I = counts(6);
N = counts(7);
O = counts(8);
P = counts(9);
S = counts(10);

X = Cl + Br + I + F;

issues = {};

% at least three atoms (water)
if sum(counts) < 3
    issues{end+1} = 'too smol';
end

% Senior Rule
senior = 2*C + N +2 - H + X;
if senior<0
    issues{end+1} = 'Fails senior rule';
end

% rings-plus-double-bonds equivalent
RDBE = C-(H+X)/2 + (N+P)/2 + 1;
if RDBE > 40
    issues{end+1} = 'Extreme RDBE';
elseif RDBE < 0
    issues{end+1} = 'negative RDBE';
end 

% element ratio checks
ratio = H/C;
if ratio <= 0.1 || ratio > 6
    issues{end+1} = 'Extreme H/C ratio';
end
ratio = F/C;
if ratio >= 6
    issues{end+1} = 'Extreme F/C ratio';
end
ratio = Cl/C;
if ratio >= 2
    issues{end+1} = 'Extreme Cl/C ratio';
end
ratio = Br/C;
if ratio >= 2
    issues{end+1} = 'Extreme Br/C ratio';
end
ratio = N/C;
if ratio >= 4
    issues{end+1} = 'Extreme N/C ratio';
end
ratio = O/C;
if ratio >= 3
    issues{end+1} = 'Extreme O/C ratio';
end
ratio = P/C;
if ratio >= 2
    issues{end+1} = 'Extreme P/C ratio';
end
ratio = S/C;
if ratio >= 3
    issues{end+1} = 'Extreme S/C ratio';
end
ratio = O/P;
if ratio < 3 % only phosphates in metabolites 
    issues{end+1} = 'non Phosphate formula';
end

%% element probability check
% NOPS
if N > 1 && O > 1 && P > 1 && S > 1
    if N < 10 || O < 20 || P < 4 || S < 3
        issues{end+1} = 'unlikly NOPS count';
    end
end

% NOP
if N > 3 && O > 3 && P > 3
    if N < 11 || O < 22 || P < 6
        issues{end+1} = 'unlikly NOP count';
    end
end
% OPS
if O > 1 && P > 1 && S > 1
    if O < 14 || P < 3 || S < 3
        issues{end+1} = 'unlikly OPS count';
    end
end
% PSN
if P > 1 && S > 1 && N > 1
    if P < 3 || S < 3 || N < 4
        issues{end+1} = 'unlikly PSN count';
    end
end
% NOS
if N > 6 && O > 6 && S > 6
    if N < 19 || O < 14 || S < 8
        issues{end+1} = 'unlikly OSN count';
    end
end

% Valid?
isValid = isempty(issues);
end
