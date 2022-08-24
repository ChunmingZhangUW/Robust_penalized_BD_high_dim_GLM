
function [hat_beta, converge] = robust_GLM_BD_parameter_estimate(...
    I_loss, family, link, X, y, beta_initial, ...
    options, index_robust_c, choice_rho_function, c_tune_constant, ...
    w_robust)

%--------------------------------------------------------------------------
% Name     : robust_GLM_BD_parameter_estimate.m
% Function : find robust version of minimum Bregman Divergence
%            parametric estimator of beta
% Criterion: 1/n*\sum_{i=1}^n \rho_q(Y_i, F^{-1}(X_i^T \bbeta))
%            w(X_i) for beta.
% Model    : F(m(x)) = beta_1 x_1 + ... + beta_K x_K,
%            where m(x)=E(Y | X=x)
% Loss     : deviance loss, exponential loss, quadratic loss,
%            quasi-likelihood
% Link F   : the canonical link function
% Algorithm: weighted quadratic approximation;
%            direct soltion, coordinate descent algorithm
% Called   : GLM_BD_parameter_estimate.m,
%            robust_weight_response_BD.m, robust_p_1_p_2_BD.m
%--------------------------------------------------------------------------
% <Input>
%    I_loss   : choice of loss function:
%               1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%               211 (V(x)=phi x)
%    family   : 0 (Gaussian), 1 (Bernoulli), 21 (Poisson_quasi)
%     link    : type of link function, 'iden', 'logit', 'log'
%      X      : n\times K design matrix, where K <= n
%      y      : n\times 1 response vector
% beta_initial: K\times 1, initial value of beta
%    options  : set of parameters
% index_robust_c:  0: \psi(r)  =  r;
%                  1: \psi(r) \ne r;
%choice_rho_function: 1 Huber \rho function; 2 Tukey biweight function
% c_tune_constant: constant used in \psi(r)
%--------------------------------------------------------------------------
% <Output>
%   hat_beta  : final estimate of beta
%   converge  : indicator of true or false convergence of the iterations
%--------------------------------------------------------------------------

if (index_robust_c == 0 || c_tune_constant == inf) && var(w_robust) == 0
    % \psi(r) = r, and w(X_i) = 1

    [hat_beta, converge] = GLM_BD_parameter_estimate(...
        I_loss, family, link, X, y, beta_initial, options);
    return
end

% (index_robust_c == 1 && c_tune_constant ~= inf) || w(x)\ne 1
% \psi(r) \ne r, or w(x)\ne 1

%%%%%%%%%%%%%%%%%%%%%%%%% for K  > n_obs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_obs, K] = size(X);

if K > n_obs; disp('K > n_obs, fails'); return; end

%%%%%%%%%%%%%%%%%%%%%%%%% for K <= n_obs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if index_robust_c == 0 && I_loss == 3 && strcmpi(link, 'iden')
    % \psi(r) = r, quadratic loss, identity link;
    % 1/n*\sum_{i=1}^n (Y_i - \X_i^T \bbeta)^2 w(X_i)

    increase = options.eps*eye(K);

    WX = zeros(n_obs, K);
    for j = 1:K
        WX(:, j) = w_robust.*X(:, j);
    end
    hat_beta = pinv(WX'*X + increase)*(WX'*y); % K*1 vector

    converge = true;
    return;
end

% index_robust_c ~= 0 || I_loss ~= 3 || strcmpi(link, 'iden')=false

%=================== step-1 initial estimate of beta ========

hat_beta = beta_initial;

%-----------------------------------------------

converge = false;
iter = 0;

if     options.robust_GLM_BD_parameter_estimate ==  1
    increase = options.eps*eye(K);

elseif options.robust_GLM_BD_parameter_estimate == 22 || ...
        options.robust_GLM_BD_parameter_estimate == 23
    square_X = X.^2;
end
while converge == false && iter <= options.maxit && ...
        max(abs(hat_beta)) <= options.max_abs_hat_beta
    iter = iter + 1;

    beta_old = hat_beta; % based on X

    %------------ weighted quadratic approximation ----------

    Xhat_beta = X*hat_beta;

    if     options.robust_GLM_BD_parameter_estimate ==  1 || ...
            options.robust_GLM_BD_parameter_estimate == 22
        [p_2_vector, p_1_d_p_2] = robust_weight_response_BD(...
            I_loss, family, link, Xhat_beta, y, options, ...
            index_robust_c, choice_rho_function, c_tune_constant);
        % p_2_vector: has been set to be >= options.delta

        ss = p_2_vector;
        zz = Xhat_beta - p_1_d_p_2;

        weight_loss = (ss.*w_robust)/n_obs;  % n_obs*1

    elseif options.robust_GLM_BD_parameter_estimate == 21
        % coordinate descent algorithm, with p_1_vector and p_2_vector

        [p_1_vector, p_2_vector] = robust_p_1_p_2_BD(...
            I_loss, family, link, Xhat_beta, y, options, ...
            index_robust_c, choice_rho_function, c_tune_constant);
        % p_2_vector: has been set to be >= options.delta

        %[p_1_vector, p_2_vector] = robust_p_1_p_2_BD_for_hat_V_n(...
        %    I_loss, family, link, hat_theta, y, options, ...
        %    index_robust_c, choice_rho_function, c_tune_constant);  % worse

        p_1_w = p_1_vector.*w_robust;
        p_2_w = p_2_vector.*w_robust;
    end

    %================= step 2: estimate of beta =============

    if     options.robust_GLM_BD_parameter_estimate ==  1
        % WQA; direct solution

        WX = zeros(n_obs, K);
        for j = 1:K
            WX(:, j) = weight_loss.*X(:, j);
        end

        beta_new = pinv(WX'*X + increase)*(WX'*zz); % K*1 vector

    elseif options.robust_GLM_BD_parameter_estimate == 21
        % WQA; coordinate descent algorithm, with p_1_vector and p_2_vector

        WX = zeros(n_obs, K);
        for j = 1:K
            WX(:, j) = p_2_w.*X(:, j);
        end
        I_2_matrix = WX'*X/n_obs;     % K*K matrix

        I_1_vector = p_1_w'*X/n_obs;  % 1*K row vector

        I_2_vector = options.eps + diag(I_2_matrix)';
        %p_2_w'*square_X/n_obs;  % 1*K row vector

        %--------------------------------------------------

        change = true;
        RUN = 0;

        beta_new = hat_beta;
        while change == true && ...
                RUN <= options.maxit && ...
                max(abs(beta_new)) <= options.max_abs_hat_beta
            RUN = RUN + 1;

            for j = 1:K
                beta_diff = (I_2_matrix(j, :)*(hat_beta - beta_new) - ...
                    I_1_vector(j))/I_2_vector(j);

                beta_new(j) = beta_diff + beta_new(j);
            end

            change = (norm(beta_new - hat_beta) > options.thresh_1);
        end

    elseif options.robust_GLM_BD_parameter_estimate == 22
        % WQA; coordinate descent algorithm, with p_2_vector and p_1_d_p_2

        wnorm2_x = options.eps + ...
            weight_loss'*square_X;  % 1*K row vector

        hat_r = zz - Xhat_beta; % residual vector

        %--------------------------------------------------

        change = true;
        RUN = 0;

        beta_new = hat_beta;
        while change == true && ...
                RUN <= options.maxit && ...
                max(abs(beta_new)) <= options.max_abs_hat_beta
            RUN = RUN + 1;

            for j = 1:K
                X_j = X(:, j);

                beta_diff = sum(weight_loss.*hat_r.*X_j)/wnorm2_x(j);

                beta_new(j) = beta_diff + beta_new(j);

                hat_r = hat_r - X_j*beta_diff;
            end

            change = (norm(beta_new - hat_beta) > options.thresh_1);
        end

    elseif options.robust_GLM_BD_parameter_estimate == 23
        % coordinate descent algorithm, with q_1_vector and q_2_vector

        change = true;
        RUN = 0;

        beta_new = hat_beta;
        while change == true && ...
                RUN <= options.maxit && ...
                max(abs(beta_new)) <= options.max_abs_hat_beta
            RUN = RUN + 1;

            for j = 1:K
                X_j        = X(:, j);
                square_X_j = square_X(:, j);

                [p_1_vector, p_2_vector] = robust_p_1_p_2_BD(...
                    I_loss, family, link, X*beta_new, y, options, ...
                    index_robust_c, choice_rho_function, c_tune_constant);
                % p_2_vector: has been set to be >= options.delta

                p_1_w = p_1_vector.*w_robust;
                p_2_w = p_2_vector.*w_robust;

                I_1_vector_j = mean(p_1_w.*X_j);
                I_2_vector_j = options.eps + mean(p_2_w.*square_X_j);
                %q_2_vector'*square_X/n_obs;  % 1*K row vector

                %--------------------------------------------------

                beta_diff = -I_1_vector_j/I_2_vector_j;

                beta_new(j) = beta_diff + beta_new(j);
            end

            change = (norm(beta_new - hat_beta) > options.thresh_1);
        end
    end

    hat_beta = beta_new;

    %--------------------- stopping rule ----------------

    converge = (norm(beta_new - beta_old) <= options.thresh_1);
end

if any(isempty(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!robust_GLM_BD_parameter_estimate.m: some estimate of beta = [ ]!!!');
end

if any(isnan(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!robust_GLM_BD_parameter_estimate.m: some estimate of beta = NaN!!!');
end

if any(isinf(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!robust_GLM_BD_parameter_estimate.m: some estimate of beta = Inf!!!');
end
