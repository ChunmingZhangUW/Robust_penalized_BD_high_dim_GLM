
function [hat_beta, converge] = GLM_BD_CD_penalized_parameter_estimate(...
    I_loss, family, link, ...
    weight_pen, penalty_set, lambda, X, y, beta_0, options)

%--------------------------------------------------------------------------
% Name     : GLM_BD_CD_penalized_parameter_estimate.m
% Function : find minimum penalized Bregman Divergence parametric estimator
% Criterion: 1/n*\sum_{i=1}^n Q(Y_i, m(X_i)) +
%            lambda \sum_{j in penalty_set} t_j |beta_j|
% Model    : F(m(x))=beta_1 x_1+...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss     : deviance loss, exponential loss, quadratic loss,
%            quasi-likelihood
% Link F   : the canonical link function
% Penalty  : L1 and weighted L1
% Algorithm: optimization: coordinate descent algorithm
% Called   : GLM_BD_parameter_estimate.m,
%            GLM_BD_CD_penalized_parameter_estimate_LM.m,
%            weight_response_BD.m, q_1_q_2_BD.m, soft_thres.m
%--------------------------------------------------------------------------
% <Input>
% I_loss : choice of loss function:
%          1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%          211 (V(x)=phi x)
% family : 0 (Gaussian); 1 (Bernoulli), 12 (Binomial); 21 (Poisson_quasi)
%  link  : type of link function, 'iden', 'logit', 'log'
%--------------------------------------------------------------------------

if lambda == 0 || isempty(penalty_set) == 1 || sum(weight_pen) == 0  % no penalty
    [hat_beta, converge] = GLM_BD_parameter_estimate(...
        I_loss, family, link, X, y, beta_0, options);
    return;
end

if I_loss == 3 && strcmpi(link, 'iden') % quadratic loss, identity link
    [hat_beta, converge] = GLM_BD_CD_penalized_parameter_estimate_LM(...
        weight_pen, penalty_set, lambda, X, y, beta_0, options);
    return;
end

[n_obs, K] = size(X);
mean_X = mean(X, 1);                         % row vector
std_X  = std(X, 1);                          % row vector

%=================== step-1 initial estimate of beta ========

hat_beta = beta_0; % based on X

%-----------------------------------------------

w = zeros(K, 1);   % K*1 column vector
w(penalty_set) = weight_pen;

Lambda = lambda*w; % K*1 column vector

zero_beta_set_1 = intersect(find(std_X == 0), find(mean_X ~= 1));
% constant columns not of the intercept
zero_beta_set_2 = find(Lambda == Inf);
% indices of j such that hat_beta(j)=0
zero_beta_set = union(zero_beta_set_1, zero_beta_set_2);
non_zero_beta_set = setdiff((1 : K)', zero_beta_set);
K_1 = length(non_zero_beta_set);

X_reduce = X(:, non_zero_beta_set);    % n_obs*K_1
Lambda_reduce = Lambda(non_zero_beta_set);

%-----------------------------------------------

converge = false;
iter = 0;

if options.GLM_BD_CD_penalized_parameter_estimate == 22
    square_X_reduce = X_reduce.^2;
end
while converge == false && iter <= options.maxit && ...
        max(abs(hat_beta)) <= options.max_abs_hat_beta
    iter = iter + 1;

    beta_old = hat_beta(non_zero_beta_set); % based on X

    %----------------------------------------------------------
    hat_theta = X_reduce*beta_old;

    if     options.GLM_BD_CD_penalized_parameter_estimate == 22
        [q_2_vector, q_1_d_q_2] = weight_response_BD(I_loss, family, link, ...
            hat_theta, y, options);
        ss = q_2_vector; % has been set to be >= options.delta
        zz = hat_theta - q_1_d_q_2;

        weight_loss = ss/n_obs;  % n_obs*1

    elseif options.GLM_BD_CD_penalized_parameter_estimate == 21
        [q_1_vector, q_2_vector] = q_1_q_2_BD(...
            I_loss, family, link, hat_theta, y, options);
        %q_2_vector has been set to be >= options.delta
    end

    %========= step 2: coordinate descent estimate of beta ==========

    if     options.GLM_BD_CD_penalized_parameter_estimate == 21
        % WQA; coordinate descent algorithm, with q_1_vector and q_2_vector

        WX = zeros(n_obs, K_1);
        for j = 1:K_1
            WX(:, j) = q_2_vector.*X_reduce(:, j);
        end
        I_2_matrix = WX'*X_reduce/n_obs;          % K_1*K_1 matrix

        I_1_vector = q_1_vector'*X_reduce/n_obs;  % 1*K_1 row vector

        I_2_vector = options.eps + diag(I_2_matrix)';
        %q_2_vector'*square_X/n_obs;  % 1*K_1 row vector

        %--------------------------------------------------

        change = true;
        RUN = 0;

        beta_new = beta_old;
        non_zero_set = 1:K_1;
        while change == true && ...
                RUN <= options.maxit && ...
                max(abs(beta_new)) <= options.max_abs_hat_beta
            RUN = RUN + 1;

            active_set = [];
            for j = non_zero_set
                beta_new(j) = soft_thres( ...
                    I_2_matrix(j, :)*(beta_old - beta_new) ...
                    - I_1_vector(j) + I_2_vector(j)*beta_new(j), ...
                    Lambda_reduce(j))/I_2_vector(j);

                if abs(beta_new(j)) > options.thresh_2
                    active_set = [active_set, j];
                else
                    beta_new(j) = 0;
                end
            end

            if options.CD_active_set == true
                non_zero_set = active_set;
            end

            change = (norm(beta_new - beta_old) > options.thresh_1);
        end

    elseif options.GLM_BD_CD_penalized_parameter_estimate == 22
        % WQA; coordinate descent algorithm, with q_2_vector and q_1_d_q_2

        wnorm2_x = options.eps + ...
            weight_loss'*square_X_reduce;  % row vector

        hat_r = zz - hat_theta; % residual vector

        %--------------------------------------------------

        change = true;
        RUN = 0;

        beta_new = beta_old;
        non_zero_set = 1:K_1;
        while change == true && ...
                RUN <= options.maxit && ...
                max(abs(beta_new)) <= options.max_abs_hat_beta
            RUN = RUN + 1;

            active_set = [];
            for j = non_zero_set
                X_reduce_j = X_reduce(:, j);

                b_old_j = beta_new(j);
                beta_new(j) = soft_thres(sum(weight_loss.*hat_r.*X_reduce_j) + ...
                    beta_new(j)*wnorm2_x(j), ...
                    Lambda_reduce(j))/wnorm2_x(j);

                if abs(beta_new(j)) > options.thresh_2
                    active_set = [active_set, j];
                else
                    beta_new(j) = 0;
                end

                beta_diff = beta_new(j) - b_old_j;

                hat_r = hat_r - X_reduce_j*beta_diff;
            end

            if options.CD_active_set == true
                non_zero_set = active_set;
            end

            change = (norm(beta_new - beta_old) > options.thresh_1);
        end
    end

    %========= step 3: final estimate of beta ==========

    hat_beta = zeros(K, 1);
    hat_beta(non_zero_beta_set) = beta_new;

    %--------------------- stopping rule ----------------

    converge = (norm(beta_new - beta_old) <= options.thresh_1);
end

if any(isempty(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!GLM_BD_CD_penalized_parameter_estimate.m: some estimate of beta = [ ]!!!');
end

if any(isnan(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!GLM_BD_CD_penalized_parameter_estimate.m: some estimate of beta = NaN!!!');
end

if any(isinf(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!GLM_BD_CD_penalized_parameter_estimate.m: some estimate of beta = Inf!!!');
end

