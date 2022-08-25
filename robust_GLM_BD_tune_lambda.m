
function [lambda_opt, min_TE, Hat_beta, I] = robust_GLM_BD_tune_lambda(...
    I_loss, family, link, X, y, beta_0, X_test, y_test, weight_pen, ...
    penalty_set, n_grid_lambda, grid_A, grid_B, misclass, options, ...
    index_robust_c, choice_rho_function, c_tune_constant, w_robust)

%--------------------------------------------------------------------------
% Name     : robust_GLM_BD_tune_lambda.m
% Function : find the optimal lambda from
%            a given training set and a given test set
% Criterion: 1/n*\sum_{i=1}^n \rho_q(Y_i, m(X_i)) w(X_i) +
%            \sum_{j in penalty_set} t_j |beta_j|
% Model    : F(m(x))=beta_1 x_1+...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss     : deviance loss, exponential loss, quadratic loss,
%            quasi-likelihood
% Link F   : the canonical link function
% Penalty  : L1 and weighted L1
% Data     : a given training set (X, y) and a given test set (X_test, y_test)
% Called   : robust_GLM_BD_CD_penalized_parameter_estimate.m,
%            robust_true_pp.m
%--------------------------------------------------------------------------
% <Input>
%  I_loss     : choice of loss function:
%               1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%               211 (V(x)=phi x)
%  family     : 0 (Gaussian), 1 (Bernoulli), 21 (Poisson_Quasi)
%   link      : type of link function, 'iden', 'logit', 'log'
%  (X, y)     : a given training set
%  beta_0     : K\times 1, initial estimate of beta for (X, y)
%  (X_test, y_test): a given test set
%  weight_pen : K_0*1 weight vector for penalized parameters, could be infinity
%  penalty_set: indices j = 1, 2, ..., K such that |beta_j| is penalized
%  n_grid_lambda     : # grid points for searching lambda
%  grid_A     : initial grid point in the log-scale
%  grid_B     : end     grid point in the log-scale
% misclass    : 1 (compute misclassification rate); 0 (compute other deviance)
%               for Bernoulli responses
%  options    : set of parameters
% index_robust_c:  0: \psi(r)  =  r;
%                  1: \psi(r) \ne r;
%choice_rho_function: 1 Huber \rho function; 2 Tukey biweight function
% c_tune_constant: constant used in \psi(r)
%     w_robust   : weight function w(X_i)
%--------------------------------------------------------------------------

K = size(X, 2);

if     options.linear_log_grid == 1
    L_grid = linspace(grid_A, grid_B, n_grid_lambda)';

elseif options.linear_log_grid == 2
    L_grid = logspace(grid_A, grid_B, n_grid_lambda)';
end

%------------------- Compute test error -------------------------
TE = zeros(n_grid_lambda, 1); Hat_beta = zeros(K, n_grid_lambda);

beta_1 = beta_0;
for k = 1:n_grid_lambda
    lambda = L_grid(k,1);

    if     options.choice_initial_value == 0
        if     options.CD_LARS == 1
            beta_1 = robust_GLM_BD_CD_penalized_parameter_estimate(...
                I_loss, family, link, weight_pen, penalty_set, lambda, X, y, ...
                beta_0, options, index_robust_c, choice_rho_function, ...
                c_tune_constant, w_robust);
        end

    elseif options.choice_initial_value == 1
        if     options.CD_LARS == 1
            beta_1 = robust_GLM_BD_CD_penalized_parameter_estimate(...
                I_loss, family, link, weight_pen, penalty_set, lambda, X, y, ...
                beta_1, options, index_robust_c, choice_rho_function, ...
                c_tune_constant, w_robust);
        end
    end
    % parameters computed from (X, y)
    Hat_beta(:, k) = beta_1;

    hat_theta_test = X_test*beta_1;
    if misclass == 1 % for Bernoulli responses
        predict_y_test = (hat_theta_test > 0);
        TE(k,1) = mean(predict_y_test ~= y_test);   % misclassification rate

    else
        predict_BD = robust_true_pp(I_loss, family, link, hat_theta_test, ...
            y_test, 0, options, index_robust_c, choice_rho_function, ...
            c_tune_constant);
        TE(k,1) = mean(predict_BD);
    end
end

[min_TE, I] = min(TE);
lambda_opt = L_grid(I);     % optimal lambda

% Hat_beta(:, I); the estimate of beta using the optimal lambda

