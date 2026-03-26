function [arrayOut, vecOut] = depadArrays(arrayIn, vecIn)
%% depadArrays Remove zero-padding rows based on a column vector.
%  depadArrays(arrayIn, vecIn) removes rows of arrayIn where the
%  corresponding entry of vecIn is zero, returning depadded versions
%  of both arrayIn and vecIn.
%
%  Inputs:
%    arrayIn : double matrix
%    vecIn   : double column vector, same number of rows as arrayIn.
%
%  Outputs:
%    arrayOut : double matrix without zero-padded rows.
%    vecOut   : double column vector without zero elements.
%
%  [arrayOut, vecOut] = depadArrays(arrayIn, vecIn)

    arguments
        arrayIn (:,:) double
        vecIn   (:,1) double {mustBeSameHeight(arrayIn, vecIn)}
    end

    % Logical index for non-zero entries in vecIn.
    idx      = (vecIn ~= 0);
    arrayOut = arrayIn(idx, :);
    vecOut   = vecIn(idx);

end  % function depadArrays


%% Local validation function
function mustBeSameHeight(matrix, vector)
% mustBeSameHeight Validate that two arrays have the same number of rows.

    if size(matrix, 1) ~= size(vector, 1)
        eid = 'Size:notEqual';
        msg = 'Inputs must have the same height (number of rows).';
        error(eid, msg)
    end

end  % function mustBeSameHeight