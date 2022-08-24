
function [hat_beta, converge] = GLM_BD_CD_penalized_parameter_estimate_SCAD(...
    I_loss, family, link, penalty_set, lambda, X, y, beta_0, options)
%--------------------------------------------------------------------------
% Name     : GLM_BD_CD_penalized_parameter_estimate_SCAD.m
% Function : find minimum penalized Bregman Divergence parametric estimator
% Criterion: 1/n*\sum_{i=1}^n Q(Y_i, m(X_i)) +
%            \sum_{j in penalty_set} P_\lambda(|beta_j|)
% Model    : F(m(x))=beta_1 x_1+...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss     : deviance loss, exponential loss; quadratic loss,
%            quasi-likelihood
% Link F   : the canonical link function
% Penalty  : SCAD
% Algorithm: LLA approximation: coordinate descent algorithm
% Called   : GLM_BD_parameter_estimate.m,
%            GLM_BD_CD_penalized_parameter_estimate_SCAD_LM.m,
%            GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA.m
%--------------------------------------------------------------------------

if lambda == 0 || isempty(penalty_set) == 1    % no penalty
    [hat_beta, converge] = GLM_BD_parameter_estimate(...
        I_loss, family, link, X, y, beta_0, options);
    return;
end

if I_loss == 3 && strcmpi(link, 'iden') % quadratic loss, identity link
    [hat_beta, converge] = GLM_BD_CD_penalized_parameter_estimate_SCAD_LM(...
        penalty_set, lambda, X, y, beta_0, options);
    return;
end

if options.CD_SCAD_LLA == 1
    [hat_beta, converge] = GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA(...
        I_loss, family, link, penalty_set, lambda, X, y, beta_0, options);
    return;
end

