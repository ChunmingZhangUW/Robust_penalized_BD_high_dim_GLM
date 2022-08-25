
function f = Poisson_pmf(J_vec, mu_vec)
%--------------------------------------------------------------------------
% Name    : Poisson_pmf.m
% Function: compute the pmf function of the Poisson distribution
% Called  : poisspdf.m
%--------------------------------------------------------------------------

f = poisspdf(J_vec, mu_vec);

return
