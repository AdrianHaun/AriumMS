function [arrayOut,vecOut] = depadArrays(arrayIn,vecIn)
%% depadArrays removes zero padding from input matrix and vector
% 
% Uses zeros in vecIn to define the padding regions and removes them from
% both inputs
%
% inputs:   arrayIn: double matrix
%           vecIn: double column vector

% outputs:  arrayOut: double matrix without padding
%           vecOut: double column vector without padding

arguments
    arrayIn (:,:) double
    vecIn   (:,1) double {mustBeSameHeight(arrayIn,vecIn)}
end

idx = vecIn ~= 0;
arrayOut = arrayIn(idx,:);
vecOut = vecIn(idx);

%% Validation function
function mustBeSameHeight(matrix,vector)
% Test for equal height
    if ~isequal(height(matrix),height(vector))
        eid = 'Size:notEqual';
        msg = 'Height of first input must equal size of second input.';
        error(eid,msg)
    end
