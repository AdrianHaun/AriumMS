function [arrayOut,vecOut] = depadArrays(arrayIn,vecIn)
%% depadArrays
% Removes zero padding from input matrix and vector

arguments
    arrayIn (:,:) double
    vecIn   (:,1) double {mustBeSameHeight(arrayIn,vecIn)}
end

idx = vecIn~=0;
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
