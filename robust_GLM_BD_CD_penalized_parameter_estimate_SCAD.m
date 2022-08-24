
function [hat_beta, converge] = ...
    robust_GLM_BD_CD_penalized_parameter_estimate_SCAD(...
    I_loss, family, link, penalty_set, lambda, X, y, beta_0, options, ...
    index_robust_c, choice_rho_function, c_tune_constant, w_robust)

%--------------------------------------------------------------------------
% Name     : robust_GLM_BD_CD_penalized_parameter_estimate_SCAD.m
% Function : find robust version of minimum penalized Bregman Divergence
%            parametric estimator
% Criterion: 1/n*\sum_{i=1}^n \rho_q(Y_i, m(X_i)) w(X_i) +
%            \sum_{j in penalty_set} P_\lambda(|beta_j|)
% Model    : F(m(x))=beta_1 x_1+...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss     : deviance loss, exponential loss; quadratic loss,
%            quasi-likelihood
% Link F   : the canonical link function
% Penalty  : SCAD
% Algorithm: LLA approximation: coordinate descent algorithm
% Called   : GLM_BD_CD_penalized_parameter_estimate_SCAD.m,
%            robust_GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA.m
%--------------------------------------------------------------------------
% <Input>
%choice_rho_function: 1 Huber \rho function; 2 Tukey biweight function
% index_robust_c:  0: \psi(r)  =  r;
%                  1: \psi(r) \ne r;
% c_tune_constant: constant used in \psi(r)
%     w_robust   : weight function w(X_i)
%--------------------------------------------------------------------------
% <Output>
%--------------------------------------------------------------------------

if (index_robust_c == 0 || c_tune_constant == inf) && var(w_robust) == 0
    % \psi(r) = r, and w(X_i) = 1

    [hat_beta, converge] = GLM_BD_CD_penalized_parameter_estimate_SCAD(...
        I_loss, family, link, penalty_set, lambda, X, y, beta_0, options);
    return
end

if options.CD_SCAD_LLA == 1
    [hat_beta, converge] = robust_GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA(...
        I_loss, family, link, penalty_set, lambda, X, y, beta_0, options, ...
        index_robust_c, choice_rho_function, c_tune_constant, w_robust);
    return
end
