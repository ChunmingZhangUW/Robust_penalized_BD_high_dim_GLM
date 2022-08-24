
function [q_1_vector, q_2_vector] = q_1_q_2_BD_for_hat_V_n(...
    I_loss, family, link, theta, y, options)

%---------------------------------------------------------------------------
% Name    : q_1_q_2_BD_for_hat_V_n.m
% Function: compute qq_1(y;\theta) and qq_2(y;\theta)
%           used in weighted quadratic approximation
%           to Bregman Divergence Q_q(y, F^{-1}(\theta))
% Model   : F(m(x))=beta_1 x_1+...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss    : deviance loss, exponential loss, quadratic loss,
%           quasi-likelihood
% Link F  : the canonical link function
% Used in : estimate the covariance matrix of the minimizer of
%           1/n*\sum_{i=1}^n Q(Y_i, m(X_i))
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

if     family == 0 && strcmpi(link, 'iden')
    % Gaussian responses, identity link

    mu = theta;

    if (I_loss == 1 || I_loss == 3) % deviance, quadratic
        q_1_vector = -2*(y - mu);

        q_2_vector = 2; %*ones(length(y), 1);
    end

elseif family == 1 && strcmpi(link, 'logit')
    % Bernoulli responses, logit link

    mu = 1./(1 + exp(-theta));

    if     I_loss == 1 % deviance
        q_1_vector = -2*(y - mu);

        if     options.UB_Bernoulli_dev_choice == 0
            q_2_vector = 2*mu.*(1 - mu);

        elseif options.UB_Bernoulli_dev_choice == 1
            q_2_vector = ones(length(y), 1)/2;
        end

    elseif I_loss == 2 % exponential
        mu_modified           = max(    mu, options.zero_thres);
        one_minus_mu_modified = max(1 - mu, options.zero_thres);

        ratio = ...
            + sqrt((1 - mu)./mu_modified).*y ...
            - sqrt(mu./one_minus_mu_modified).*(1 - y);
        q_1_vector = -ratio/2;

        q_2_vector = ( ...
            + sqrt(mu./one_minus_mu_modified).*(1 - y) ...
            + sqrt((1 - mu)./mu_modified).*y )/4;
    end

elseif family == 12 && strcmpi(link, 'logit')
    % Binomial responses, logit link

    y_Bin = y(:, 1);
    N_Bin = y(:, 2);

    p_Bin = 1./(1 + exp(-theta));

    mu = N_Bin.*p_Bin;

    V_mu_original = mu.*(1 - p_Bin);    % mu.*(1-p);

    if I_loss == 1 % deviance
        q_1_vector = -2*(y_Bin - mu);

        if     options.UB_Binomial_dev_choice == 0
            q_2_vector = 2*V_mu_original;

        elseif options.UB_Binomial_dev_choice == 1
            q_2_vector = N_Bin/2;
        end
    end

elseif family == 21 && strcmpi(link, 'log')
    % Poisson_quasi, log link

    theta(theta < -options.BD_C) = -options.BD_C;
    theta(theta >  options.BD_C) =  options.BD_C;

    mu = exp(theta);

    if I_loss == 211  % V(x)=phi x, (negative) quasi-likelihood
        q_1_vector = -(y - mu);

        q_2_vector = mu;
    end
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if any(isempty(q_1_vector)) == 1
    disp(' !!!q_1_q_2_BD_for_hat_V_n.m: some estimate of q_1_vector = [ ]!!!');
end

if any(isnan(q_1_vector))   == 1
    disp(' !!!q_1_q_2_BD_for_hat_V_n.m: some estimate of q_1_vector = NaN!!!');
end

if any(isinf(q_1_vector))   == 1
    disp(' !!!q_1_q_2_BD_for_hat_V_n.m: some estimate of q_1_vector = Inf!!!');
end
%-------------

if any(isempty(q_2_vector)) == 1
    disp(' !!!q_1_q_2_BD_for_hat_V_n.m: some estimate of q_2_vector = [ ]!!!');
end

if any(isnan(q_2_vector))   == 1
    disp(' !!!q_1_q_2_BD_for_hat_V_n.m: some estimate of q_2_vector = NaN!!!');
end

if any(isinf(q_2_vector))   == 1
    disp(' !!!q_1_q_2_BD_for_hat_V_n.m: some estimate of q_2_vector = Inf!!!');
end
