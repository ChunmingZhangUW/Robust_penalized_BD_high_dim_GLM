
function [hat_beta, converge] = ...
    robust_GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA(...
    I_loss, family, link, penalty_set, lambda, X, y, beta_0, options, ...
    index_robust_c, choice_rho_function, c_tune_constant, w_robust)

%--------------------------------------------------------------------------
% Name     : robust_GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA.m
% Function : find minimum penalized Bregman Divergence parametric estimator
% Criterion: 1/n*\sum_{i=1}^n \rho_q(Y_i, m(X_i)) w(X_i) +
%            \sum_{j in penalty_set} P_\lambda(|beta_j|)
% Model    : F(m(x))=beta_1 x_1+...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss     : deviance loss, exponential loss; quadratic loss,
%            quasi-likelihood
% Link F   : the canonical link function
% Penalty  : SCAD
% Algorithm: LLA approximation: coordinate descent algorithm
% Called   : robust_GLM_BD_parameter_estimate.m,
%            GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA.m,
%            true_penalty.m,
%            robust_weight_response_BD.m,
%            robust_p_1_p_2_BD.m, soft_thres.m,
%--------------------------------------------------------------------------
% <Input>
% index_robust_c:  0: \psi(r)  =  r;
%                  1: \psi(r) \ne r;
%choice_rho_function: 1 Huber \rho function; 2 Tukey biweight function
% c_tune_constant: constant used in \psi(r)
%     w_robust   : weight function w(X_i)
%--------------------------------------------------------------------------
% <Output>
%--------------------------------------------------------------------------

if lambda == 0 || isempty(penalty_set) == 1  % no penalty
    [hat_beta, converge] = robust_GLM_BD_parameter_estimate...
        (I_loss, family, link, ...
        X, y, beta_0, options, index_robust_c, choice_rho_function, ...
        c_tune_constant, w_robust);
    return;
end

if (index_robust_c == 0 || c_tune_constant == inf) && var(w_robust) == 0
    % \psi(r) = r, and w(X_i) = 1

    [hat_beta, converge] = GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA(...
        I_loss, family, link, penalty_set, lambda, X, y, beta_0, options);
    return
end

[n_obs, K] = size(X);
mean_X = mean(X, 1);                         % row vector
std_X  = std(X, 1);                          % row vector

%----------------------------------------------------------

I_penalty = 2; % SCAD penalty

%=================== step-1 initial estimate of beta ========

hat_beta = beta_0; % based on X

%-----------------------------------------------

converge = false;
iter = 0;

zero_beta_set_1 = intersect(find(std_X == 0), find(mean_X ~= 1));
% constant columns not of the intercept
while converge == false && iter <= options.maxit && ...
        max(abs(hat_beta)) <= options.max_abs_hat_beta
    iter = iter + 1;

    %-----------------------------------------------

    w = zeros(K, 1);   % K*1 column vector
    for j = penalty_set
        if     hat_beta(j) ~= 0
            w(j) = true_penalty(hat_beta(j), I_penalty, lambda, 1, options) ...
                * sign(hat_beta(j))/lambda;

        elseif hat_beta(j) == 0
            w(j) = 1;
        end
    end

    Lambda = lambda*w; % K*1 column vector

    zero_beta_set_2 = find(Lambda == Inf);
    % indices of j such that hat_beta(j)=0
    zero_beta_set = union(zero_beta_set_1, zero_beta_set_2);
    non_zero_beta_set = setdiff((1 : K)', zero_beta_set);
    K_1 = length(non_zero_beta_set);

    X_reduce = X(:, non_zero_beta_set);
    Lambda_reduce = Lambda(non_zero_beta_set);

    beta_old = hat_beta(non_zero_beta_set); % based on X

    %-----------------------------------------------

    hat_theta = X_reduce*beta_old;

    if     options.robust_GLM_BD_CD_penalized_parameter_estimate == 22
        square_X_reduce = X_reduce.^2;
        [p_2_vector, p_1_d_p_2] = robust_weight_response_BD(...
            I_loss, family, link, hat_theta, y, options, ...
            index_robust_c, choice_rho_function, c_tune_constant);
        ss = p_2_vector.*w_robust; % has been set to be >= options.delta
        zz = hat_theta - p_1_d_p_2;

        weight_loss = ss/n_obs;  % n_obs*1

    elseif options.robust_GLM_BD_CD_penalized_parameter_estimate == 21
        [p_1_vector, p_2_vector] = robust_p_1_p_2_BD(...
            I_loss, family, link, hat_theta, y, options, ...
            index_robust_c, choice_rho_function, c_tune_constant);
        %p_2_vector has been set to be >= options.delta

        p_1_w = p_1_vector.*w_robust;
        p_2_w = p_2_vector.*w_robust;
    end

    %========= step 2: coordinate descent estimate of beta ==========

    if     options.robust_GLM_BD_CD_penalized_parameter_estimate == 21
        % WQA; coordinate descent algorithm, with p_1_vector and p_2_vector

        WX = zeros(n_obs, K);
        for j = 1:K_1
            WX(:, j) = p_2_w.*X_reduce(:, j);
        end
        I_2_matrix = WX'*X_reduce/n_obs;     % K_1*K_1 matrix

        I_1_vector = p_1_w'*X_reduce/n_obs;  % 1*K_1 row vector

        I_2_vector = options.eps + diag(I_2_matrix)';
        %p_2_w'*square_X/n_obs;  % 1*K_1 row vector

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

    elseif options.robust_GLM_BD_CD_penalized_parameter_estimate == 22
        % WQA; coordinate descent algorithm, with q_2_vector and q_1_d_q_2

        wnorm2_x = options.eps + ...
            weight_loss'*square_X_reduce; % row vector

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
    disp(' !!!robust_GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA.m: some estimate of beta = [ ]!!!');
end

if any(isnan(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!robust_GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA.m: some estimate of beta = NaN!!!');
end

if any(isinf(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!robust_GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA.m: some estimate of beta = Inf!!!');
end

