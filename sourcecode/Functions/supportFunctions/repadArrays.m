function [arrayOut, vecOut] = repadArrays(arrayIn, vecIn, padSize)
%% repadArrays Add zero padding to a matrix and a vector to match a given height.
%
% repadArrays pads the input matrix arrayIn and column vector vecIn with zeros
% at the bottom so that the number of rows of the outputs equals padSize.
%
% Inputs:
%   arrayIn  - double matrix. Number of rows defines original height.
%   vecIn    - double column vector. Must have same number of rows as arrayIn.
%   padSize  - scalar integer specifying desired number of rows for outputs.
%
% Outputs:
%   arrayOut - double matrix padded with zeros to have padSize rows.
%   vecOut   - double column vector padded with zeros to have padSize rows.
%
% The function assumes padSize is greater than or equal to the height of
% the inputs. Use repadArrays when you need to align array and vector
% sizes by zero-padding them to a common row count.

arguments
    arrayIn (:,:) double
    vecIn   (:,1) double {mustBeSameHeight(arrayIn, vecIn)}
    padSize (1,1) {mustBeInteger, mustBeBiggerThan(padSize, vecIn)}
end

originalHeight = size(arrayIn, 1);

arrayOut = padarray(arrayIn, padSize - originalHeight, 0, 'post');
vecOut   = padarray(vecIn,   padSize - originalHeight, 0, 'post');

end


%% Validation functions
function mustBeSameHeight(matrix, vector)
% mustBeSameHeight Validate that matrix and vector have the same height.
    if ~isequal(size(matrix, 1), size(vector, 1))
        eid = 'Size:notEqual';
        msg = 'Height (number of rows) of the first input must equal that of the second input.';
        error(eid, msg)
    end
end


function mustBeBiggerThan(inputValue, vector)
% mustBeBiggerThan Validate that padSize is not smaller than the vector height.
    if inputValue < size(vector, 1)
        eid = 'Value:tooLow';
        msg = 'padSize must be greater than or equal to the height of the inputs.';
        error(eid, msg)
    end
end