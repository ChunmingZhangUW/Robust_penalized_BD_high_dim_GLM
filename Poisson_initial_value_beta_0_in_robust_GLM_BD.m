
function [beta_0, converge] = ...
    Poisson_initial_value_beta_0_in_robust_GLM_BD(...
    I_loss, family, link, X, y, options)
%
% Name    : Poisson_initial_value_beta_0_in_robust_GLM_BD.m
% Function: obtain initial value of beta_0 in robust-BD GLM
% Called  : weight_function_in_robust_BD_estimation.m, 
%           robust_GLM_BD_parameter_estimate.m
%--------------------------------------------------------------

[n_obs, K] = size(X);

%=====================================================================

index_robust_o = [1 1];
choice_rho_function_o = options.choice_rho_function;
c_tune_constant_o     = options.c_tune_constant;

w_robust_o = weight_function_in_robust_BD_estimation...
    (X, options, index_robust_o(2));

%======== initial estimate of \beta_0 =========================

mean_y = mean(y);
beta_0 = [log(mean_y+0.1); zeros(K-1, 1)];   % starting value of beta

[beta_0] = robust_GLM_BD_parameter_estimate(...
    I_loss, family, link, X, y, beta_0, options, ...
    0, choice_rho_function_o, ...
    c_tune_constant_o, ones(n_obs, 1));
% iterative solution from classical-BD

% beta_0

%=============== test the initial estimates ==============================

% robust estimation
[~, converge] = robust_GLM_BD_parameter_estimate(...
    I_loss, family, link, X, y, beta_0, options, index_robust_o(1), ...
    choice_rho_function_o, c_tune_constant_o, w_robust_o);

if converge == false
    disp('Poisson_initial_value_beta_0_in_robust_GLM_BD.m: converge = false')
end