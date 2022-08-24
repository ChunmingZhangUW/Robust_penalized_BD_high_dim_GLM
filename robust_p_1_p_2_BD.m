
function [p_1_vector, p_2_vector, results] = robust_p_1_p_2_BD(...
    I_loss, family, link, theta, y, options, ...
    index_robust_c, choice_rho_function, c_tune_constant)

%--------------------------------------------------------------------------
% Name    : robust_p_1_p_2_BD.m
% Function: find pp_1(y;\theta) and pp_2(y;\theta)
%           used in weighted quadratic approximation
%           to robust Bregman Divergence \rho_q(y, F^{-1}(\theta))
% Model   : F(m(x))=beta_1 x_1+...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss    : deviance loss, exponential loss, quadratic loss,
%           quasi-likelihood
% Link F  : the canonical link function
% Used in : LARS algorithm
% Called  : q_1_q_2_BD.m, robust_p_1_p_2_BD_for_hat_V_n.m
%--------------------------------------------------------------------------
% <Input>
%    I_loss   : choice of loss function:
%               1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%               211 (V(x)=phi x)
%    family   : 0 (Gaussian); 1 (Bernoulli), 12 (Binomial); 21 (Poisson_quasi)
%     link    : type of link function, 'iden', 'logit', 'log'
%    theta    : scalar or vector
%      y      : scalar or vector with the same size of theta
%    beta_0   : K\times 1, initial value of beta
%    options  : set of parameters
% index_robust_c:  0: \psi(r)  =  r;
%                  1: \psi(r) \ne r;
%choice_rho_function: 1 Huber \rho function; 2 Tukey biweight function
% c_tune_constant: constant used in \psi(r)
%--------------------------------------------------------------------------
% <Output>
%  p_1_vector : same size of y
%  p_2_vector : same size of y
%--------------------------------------------------------------------------

%--------------------- for classical-BD ----------------------

if (index_robust_c == 0 || c_tune_constant == inf)   % \psi(r) = r

    [p_1_vector, p_2_vector] = q_1_q_2_BD(...
        I_loss, family, link, theta, y, options);
    return
end

%---------------------- for robust-BD ------------------------

% index_robust_c == 1 && c_tune_constant ~= inf  % \psi(r) \ne r

[p_1_vector, p_2_vector, results] = robust_p_1_p_2_BD_for_hat_V_n(...
    I_loss, family, link, theta, y, options, ...
    index_robust_c, choice_rho_function, c_tune_constant);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_2_vector = max(p_2_vector, options.delta);  % set to be >= options.delta

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if any(isempty(p_1_vector)) == 1
    disp(' !!!robust_p_1_p_2_BD.m: some estimate of p_1_vector = [ ]!!!');
end

if any(isnan(p_1_vector))   == 1
    disp(' !!!robust_p_1_p_2_BD.m: some estimate of p_1_vector = NaN!!!');
end

if any(isinf(p_1_vector))   == 1
    disp(' !!!robust_p_1_p_2_BD.m: some estimate of p_1_vector = Inf!!!');
end
%-------------

if any(isempty(p_2_vector)) == 1
    disp(' !!!robust_p_1_p_2_BD.m: some estimate of p_2_vector = [ ]!!!');
end

if any(isnan(p_2_vector))   == 1
    disp(' !!!robust_p_1_p_2_BD.m: some estimate of p_2_vector = NaN!!!');
end

if any(isinf(p_2_vector))   == 1
    disp(' !!!robust_p_1_p_2_BD.m: some estimate of p_2_vector = Inf!!!');
end