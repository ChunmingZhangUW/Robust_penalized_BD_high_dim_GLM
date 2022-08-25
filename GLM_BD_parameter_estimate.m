
function [hat_beta, converge, iter, RUN] = GLM_BD_parameter_estimate(...
    I_loss, family, link, X, y, beta_0, options)

%--------------------------------------------------------------------------
% Name     : GLM_BD_parameter_estimate.m
% Function : find minimum Bregman Divergence parametric estimator
% Criterion: 1/n_obs*\sum_{i=1}^n_obs Q_q(Y_i, m(X_i))
% Model    : F(m(x)) = beta_1 x_1 + ... + beta_K x_K, where m(x)=E(Y|X=x)
% Loss     : deviance loss, exponential loss, quadratic loss,
%            quasi-likelihood
% Link F   : the canonical link function
% Algorithm: weighted quadratic approximation;
%            direct soltion, coordinate descent algorithm
% Called   : weight_response_BD.m, q_1_q_2_BD.m
%--------------------------------------------------------------------------
% <Input>
%    I_loss   : choice of loss function:
%               1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%               211 (V(x)=phi x)
%    family   : 0 (Gaussian), 1 (Bernoulli), 21 (Poisson_quasi)
%     link    : type of link function, 'iden', 'logit', 'log'
%      X      : n_obs\times K design matrix, where K <= n_obs
%      y      : n_obs\times 1 response vector
%    beta_0   : K\times 1, initial value of beta
%    options  : set of parameters
%--------------------------------------------------------------------------
% <Output>
%   hat_beta  : final estimate of beta
%   converge  : indicator of true or false convergence of the iterations
%--------------------------------------------------------------------------

[n_obs, K] = size(X);

if K > n_obs; disp('K > n_obs, fails'); return; end

%%%%%%%%%%%%%%%%%%%%%%%%% for K <= n_obs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if I_loss == 3 && strcmpi(link, 'iden')
    % quadratic loss, and identity link
    increase = options.eps*eye(K);

    hat_beta = (X'*X + increase)\(X'*y); % K*1 vector

    converge = true;
    iter = 0;
    RUN = 0;
    return;
end

%=================== step-1 initial estimate of beta ========

hat_beta = beta_0;

%-----------------------------------------------

converge = false;
iter = 0;

if     options.GLM_BD_parameter_estimate ==  1
    increase = options.eps*eye(K);

elseif options.GLM_BD_parameter_estimate == 22 || ...
        options.GLM_BD_parameter_estimate == 23
    square_X = X.^2;
end
while converge == false && iter <= options.maxit && ...
        max(abs(hat_beta)) <= options.max_abs_hat_beta
    iter = iter + 1;

    beta_old = hat_beta; % based on X

    %------------ weighted quadratic approximation ----------

    hat_theta = X*hat_beta;

    if     options.GLM_BD_parameter_estimate ==  1 || ...
            options.GLM_BD_parameter_estimate == 22
        [q_2_vector, q_1_d_q_2] = weight_response_BD(...
            I_loss, family, link, hat_theta, y, options);
        % q_2_vector: has been set to be >= options.delta

        ss = q_2_vector;
        zz = hat_theta - q_1_d_q_2;

        weight_loss = ss/n_obs;  % n_obs*1

    elseif options.GLM_BD_parameter_estimate == 21
        % coordinate descent algorithm, with q_1_vector and q_2_vector

        [q_1_vector, q_2_vector] = q_1_q_2_BD(...
            I_loss, family, link, hat_theta, y, options);
        % q_2_vector: has been set to be >= options.delta

        %[q_1_vector, q_2_vector] = q_1_q_2_BD_for_hat_V_n(...
        %    I_loss, family, link, hat_theta, y, options);  % worse
    end

    %================= step 2: estimate of beta =============

    if     options.GLM_BD_parameter_estimate ==  1
        % WQA; direct solution

        WX = zeros(n_obs, K);
        for j = 1:K
            WX(:, j) = weight_loss.*X(:, j);
        end

        beta_new = pinv(WX'*X + increase)*(WX'*zz); % K*1 vector

        RUN = 0;

    elseif options.GLM_BD_parameter_estimate == 21
        % WQA; coordinate descent algorithm, with q_1_vector and q_2_vector

        WX = zeros(n_obs, K);
        for j = 1:K
            WX(:, j) = q_2_vector.*X(:, j);
        end
        I_2_matrix = WX'*X/n_obs;          % K*K matrix

        I_1_vector = q_1_vector'*X/n_obs;  % 1*K row vector

        I_2_vector = options.eps + ...
            diag(I_2_matrix)';  % 1*K row vector

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

    elseif options.GLM_BD_parameter_estimate == 22
        % WQA; coordinate descent algorithm, with q_2_vector and q_1_d_q_2

        wnorm2_x = options.eps + ...
            weight_loss'*square_X;  % 1*K row vector

        hat_r = zz - hat_theta; % residual vector

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

    elseif options.GLM_BD_parameter_estimate == 23 % 
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

                [q_1_vector, q_2_vector] = q_1_q_2_BD(...
                    I_loss, family, link, X*beta_new, y, options);
                % q_2_vector: has been set to be >= options.delta

                I_1_vector_j = mean(q_1_vector.*X_j);
                I_2_vector_j = options.eps + mean(q_2_vector.*square_X_j);
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
    disp(' !!!GLM_BD_parameter_estimate.m: some estimate of beta = [ ]!!!');
end

if any(isnan(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!GLM_BD_parameter_estimate.m: some estimate of beta = NaN!!!');
end

if any(isinf(hat_beta)) == 1
    disp(['I_loss = ', num2str(I_loss)]);
    disp(' !!!GLM_BD_parameter_estimate.m: some estimate of beta = Inf!!!');
end

%%%%%%%%%%%%%%%%%%%%%%%% Notes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab file: GLM_BD_parameter_estimate.m
%
% option =  1: direct method by Newton-Raphson iteration
%
% option = 21: CD method, using q_1_vector and q_2_vector
%
% option = 22: CD method, using q_1_vector and q_1_vector/q_2_vector
%              (original version)
%
% option = 23: CD method, using q_1_vector and q_2_vector
%