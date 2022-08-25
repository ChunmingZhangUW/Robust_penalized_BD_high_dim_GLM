
function F = Poisson_cdf(J_vec, mu_vec)
%--------------------------------------------------------------------------
% Name    : Poisson_cdf.m
% Function: compute the cdf function of the Poisson distribution
% Called  : poisscdf.m
%--------------------------------------------------------------------------

F = poisscdf(J_vec, mu_vec);

return
