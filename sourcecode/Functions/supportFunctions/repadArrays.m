function [arrayOut,vecOut] = repadArrays(arrayIn,vecIn,padSize)
%% repadArrays
% adds zero padding to array and vector to match padSize
%
% inputs:   arrayIn: double matrix
%           vecIn: double column vector
%           padSize: height of outputs
% outputs:  arrayOut: double matrix padded with zeros to height of padSize
%           vecOut: double column vector padded with zeros to height of padSize


arguments
    arrayIn (:,:) double
    vecIn   (:,1) double {mustBeSameHeight(arrayIn,vecIn)}
    padSize (1,1) {mustBeInteger,mustBeBiggerThan(padSize,vecIn)}
end

orginalHeight = height(arrayIn);
arrayOut = padarray(arrayIn,padSize-orginalHeight,0,'post');
vecOut = padarray(vecIn,padSize-orginalHeight,0,'post');

%% Validation functions
function mustBeSameHeight(matrix,vector)
% Test for equal height
    if ~isequal(height(matrix),height(vector))
        eid = 'Size:notEqual';
        msg = 'Height of first input must equal size of second input.';
        error(eid,msg)
    end

function mustBeBiggerThan(input,vector)
% Test if padSize is set correctly
    if input < height(vector)
        eid = 'Value:tooLow';
        msg = 'padSize must be greater than height of other inputs';
        error(eid,msg)
    end
