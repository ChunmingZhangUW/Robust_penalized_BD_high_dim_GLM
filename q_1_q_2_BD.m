
function [q_1_vector, q_2_vector] = q_1_q_2_BD(...
    I_loss, family, link, theta, y, options)

%--------------------------------------------------------------------------
% Name    : q_1_q_2_BD.m
% Function: compute qq_1(y;\theta) and qq_2(y;\theta)
%           used in weighted quadratic approximation
%           to Bregman Divergence Q_q(y, F^{-1}(\theta))
% Model   : F(m(x))=beta_1 x_1+...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss    : deviance loss, exponential loss, quadratic loss,
%           quasi-likelihood
% Link F  : the canonical link function
% Called  : q_1_q_2_BD_for_hat_V_n.m
%--------------------------------------------------------------------------
% <Input>
% I_loss : choice of loss function:
%          1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%          211 (V(x)=phi x)
% family : 0 (Gaussian); 1 (Bernoulli), 12 (Binomial); 21 (Poisson_quasi)
%  link  : type of link function, 'iden', 'logit', 'log'
%  theta : scalar or vector
%    y   : scalar or vector with the same size of theta
% options: set of parameters
%--------------------------------------------------------------------------
% <Output>
% q_1_vector: same size of y
% q_2_vector: same size of y
%--------------------------------------------------------------------------

[q_1_vector, q_2_vector] = q_1_q_2_BD_for_hat_V_n(...
    I_loss, family, link, theta, y, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q_2_vector = max(q_2_vector, options.delta);  % set to be >= options.delta

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if any(isempty(q_1_vector)) == 1
    disp(' !!!q_1_q_2_BD.m: some estimate of q_1_vector = [ ]!!!');
end

if any(isnan(q_1_vector))   == 1
    disp(' !!!q_1_q_2_BD.m: some estimate of q_1_vector = NaN!!!');
end

if any(isinf(q_1_vector))   == 1
    disp(' !!!q_1_q_2_BD.m: some estimate of q_1_vector = Inf!!!');
end
%-------------

if any(isempty(q_2_vector)) == 1
    disp(' !!!q_1_q_2_BD.m: some estimate of q_2_vector = [ ]!!!');
end

if any(isnan(q_2_vector))   == 1
    disp(' !!!q_1_q_2_BD.m: some estimate of q_2_vector = NaN!!!');
end

if any(isinf(q_2_vector))   == 1
    disp(' !!!q_1_q_2_BD.m: some estimate of q_2_vector = Inf!!!');
end
