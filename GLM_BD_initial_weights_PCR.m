
function [beta, iter] = GLM_BD_initial_weights_PCR(...
    I_loss, family, link, X, y, kappa, options)

%--------------------------------------------------------------------------
% Name     : GLM_BD_initial_weights_PCR.m
% Function : find initial weights via componentwise regression
% Criterion: 1/n*\sum_{i=1}^n Q(Y_i, m(X_{ij})) + kappa |\beta_j|,
% Model    : F(m(X_{ij}))=\alpha_j+X_{ij} \beta_j, where m(x)=E(Y|X=x)
% Loss     : deviance loss, exponential loss, quadratic loss,
%            quasi-likelihood
% Link F   : the canonical link function
% Penalty  : L1 and weighted L1
% Algorithm: coordinate descent
% Called   : soft_thres.m, weight_response_BD.m
%--------------------------------------------------------------------------
% <Input>
%   I_loss      : choice of loss function:
%                 1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%                 211 (V(x)=phi x)
%   family      : 0 (Gaussian), 1 (Bernoulli), 21 (Poisson_quasi)
%    link       : type of link function, 'iden', 'logit', 'log'
%     X         : n\times K design matrix
%     y         : n\times 1 response vector
%   kappa       : regularization parameter, a scalar
%   options     : set of parameters
%--------------------------------------------------------------------------
% <Output>
%    beta       : K\times 1, PCR estimates of beta_j
%--------------------------------------------------------------------------

[n_obs, K] = size(X);

mean_X = mean(X, 1);                        % row vector
std_X  = std(X, 1);                         % row vector
mean_y = mean(y);

%-----------------------------------------------

sig_set = find(std_X > 0); % indices of significant variables

beta = Inf(K, 1); % based on X, correspond to weight = 0

if I_loss == 3 && strcmpi(link, 'iden')  % quadratic loss, identity link
    delta = eps;

    center_X_sig_set = X(:, sig_set) - ones(n_obs,1)*mean_X(sig_set);
    center_y = y - mean_y;
    beta(sig_set) = soft_thres((2/n_obs)*center_X_sig_set'*center_y, kappa) ...
        ./((2/n_obs)*sum(center_X_sig_set.^2, 1) + delta)';

    %     if beta == 0
    %         disp(' !!!initial_weights_PCR.m: all estimates of beta = 0!!!');
    %     end

    if any(isnan(beta)) == 1
        disp(' !!!initial_weights_PCR.m: some estimate of beta = NaN!!!');
    end

    if any(isinf(beta)) == 1
        disp(' !!!initial_weights_PCR.m: some estimate of beta = Inf!!!');
    end

    iter = 1;
    return;
end

%-----------------------------------------------

if     family ==  0  % Gaussian responses
    a_0 = mean_y;

elseif family ==  1   % Bernoulli responses
    a_0 = log((mean_y+0.1)/(1-mean_y+0.1));

elseif family == 12  % Binomial responses
    y_Bin = y(:, 1);
    N_Bin = y(:, 2);

    a_0 = log((mean(y_Bin)+eps)/(mean(N_Bin-y_Bin)+eps));

elseif family == 21  % Poisson quasi-likelihood model
    a_0 = log(mean_y+0.1);
end

for j = sig_set
    %--------------------- centralize and standardize X -----

    x_j_star = (X(:, j) - mean_X(j))/std_X(j);  % vector
    kappa_j_star = kappa/std_X(j);

    %=================== step-1 initial estimate of beta ========

    b_old_j = 0;        % based on X
    a_old_j = a_0;      % F(E(Y))

    b_old_j_star = b_old_j*std_X(j);   % based on X^*
    a_old_j_star = a_old_j + mean_X(j)*b_old_j;

    hat_a_b_j_star = [a_old_j_star; b_old_j_star]; % 2*1 vector

    %=================== step-2 iterative estimate of beta ========

    converge = false;
    iter = 0;

    while converge == false && iter <= options.maxit && ...
            max(abs(hat_a_b_j_star)) <= options.max_abs_hat_beta
        iter = iter + 1;

        b_old_j_star = hat_a_b_j_star(2);
        a_old_j_star = hat_a_b_j_star(1);

        %------------ weighted quadratic approximation ----------

        hat_theta_j = hat_a_b_j_star(1) + x_j_star*hat_a_b_j_star(2);

        [q_2_vector, q_1_d_q_2] = weight_response_BD(I_loss, family, link, ...
            hat_theta_j, y, options);
        ss = q_2_vector; % has been set to be >= options.delta
        zz = hat_theta_j - q_1_d_q_2;

        %------------- estimate of alpha, beta ------------------

        weight_loss = ss/n_obs;  % n_obs*1
        sum_weight_loss = sum(weight_loss);

        wmean_z = sum(weight_loss.*zz)/sum_weight_loss;
        wmean_X_j_star = sum(weight_loss.*x_j_star)/sum_weight_loss;

        wcenter_X_j_star = x_j_star - wmean_X_j_star;
        wnorm2_center_X_j_star = options.delta + ...
            weight_loss'*(wcenter_X_j_star.^2);

        b_new_j_star = soft_thres(...
            sum(weight_loss.*(zz - wmean_z).*wcenter_X_j_star), ...
            kappa_j_star)/wnorm2_center_X_j_star;
        a_new_j_star = wmean_z - wmean_X_j_star*b_new_j_star;

        hat_a_b_j_star = [a_new_j_star; b_new_j_star]; % 2*1 vector

        %===================== stopping rule ==============================

        converge = (norm([a_new_j_star; b_new_j_star] - ...
            [a_old_j_star; b_old_j_star]) <= options.thresh_1);
    end

    %=================== step 3: re-scaled estimate of beta ==========

    beta(j) = b_new_j_star./std_X(j); % based on X
end

% if beta == 0
%     disp(' !!!initial_weights_PCR.m: all estimates of beta = 0!!!');
% end

if any(isnan(beta)) == 1
    disp(' !!!GLM_BD_initial_weights_PCR.m: some estimate of beta = NaN!!!');
end

if any(isinf(beta)) == 1
    disp(' !!!GLM_BD_initial_weights_PCR.m: some estimate of beta = Inf!!!');
end

