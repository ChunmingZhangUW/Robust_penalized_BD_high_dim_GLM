
function b = soft_thres(x, r)
%--------------------------------------------------------------------------
% Name    : soft_thres.m
% Function: soft thresholding operator
%--------------------------------------------------------------------------

abs_x = abs(x);

b = sign(x) .* max(abs_x - r, 0);
