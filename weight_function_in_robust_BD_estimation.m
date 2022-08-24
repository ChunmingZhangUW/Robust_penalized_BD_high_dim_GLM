
function w_robust = weight_function_in_robust_BD_estimation(...
    tilde_X, options, index_robust_w)

%--------------------------------------------------------------------------
% Name    : weight_function_in_robust_BD_estimation.m
% Function: compute weight function w(x) for robust-BD estimation
%--------------------------------------------------------------------------
% <Input>
%    tilde_X    : n_obs*K design matrix
%    options    :
% index_robust_w: choice of weight function for the X-space;
%                 0: w(x)  =  1;
%                 1: w(x) \ne 1;
%                 index_robust = [index_robust_c index_robust_w]
%--------------------------------------------------------------------------
% <Output>
%   w_robust    : weight vector, n_obs*1
%--------------------------------------------------------------------------

[n_obs, K] = size(tilde_X);
vector_one = ones(n_obs, 1); % vector (1,...,1)^T, n_obs*1

%-------------------- tilde_X is a constant column ------------------------
if K == 1 && var(tilde_X) == 0  % tilde_X is a constant column
    w_robust = vector_one;
    return
end

%------------------------------- K >= 2 -----------------------------------

X = tilde_X(:, 2:end); % excluding column of the intercept term
[p_n] = size(X, 2);

if     index_robust_w == 0    % w(x_1,...,x_{p_n})  =  1
    w_robust = vector_one;
    return

elseif index_robust_w == 1     % w(x_1,...,x_{p_n}) \ne 1
    digit = num2str(options.choice_of_weight_function_in_robust_BD_estimation);
    digit_1 = digit(1);

    if     digit_1 == '1' % for any p_n >= 1
        % w(x_{i,1},...,x_{i,p_n}) = 1/\sqrt{1+\sum_{j=1}^{p_n} ((x_{i,j}-m_j)/s_j)^2}

        %------------------- median vector --------------------------------
        M_n = median(X, 1)'; % median vector,   p_n*1, M_n = (m_1,...,m_{p_n})^T

        %----------------- variance vector --------------------------------
        S_n = zeros(p_n, 1); % variance vector, p_n*1, S_n = (s_1,...,s_{p_n})^T

        if     options.choice_of_weight_function_in_robust_BD_estimation == 10
            % S_n = 1
            S_n = ones(p_n, 1);

        elseif options.choice_of_weight_function_in_robust_BD_estimation == 11
            % MAD estimate of scale,
            % S_n(j, 1) = median_i( |x_{ij} - m_j|: i=1,...,n )
            %S_n = median( abs(X - vector_one*M_n'), 1)';  % p_n*1

            for j = 1:p_n
                X_j = X(:, j);
                S_n(j, 1) = median( abs(X_j - M_n(j)) );
            end

        elseif options.choice_of_weight_function_in_robust_BD_estimation == 12
            % Rousseeuw and Croux (1993),
            % S_n(j, 1) = 1.1926*...
            % median_i{ median_k( |x_{ij} - x_{kj}|: k=1,...,n ): i=1,...,n }
            % very slow
            med_in = zeros(n_obs, 1);

            for j = 1:p_n
                X_j = X(:, j);

                for i = 1:n_obs
                    med_in(i) = median( abs(X_j(i) - X_j) );
                end
                S_n(j, 1) = 1.1926*median( med_in );
            end
        end
        for j = 1:p_n
            S_n(j) = max(S_n(j), eps);
        end

        %---------------- standardized design matrix X --------------------
        std_X = (X - vector_one*M_n')./(vector_one*S_n');

        w_robust = 1./sqrt(1 + sum(std_X.^2, 2));

    elseif options.choice_of_weight_function_in_robust_BD_estimation == 2 ...
            && p_n <= n_obs
        % w(x_{i,1},...,x_{i,p_n}) = \sqrt{1 - H(i, i)}
        Hat_matrix = X*pinv(X'*X)*X';

        w_robust = sqrt(1 - diag(Hat_matrix));
        return

    elseif digit_1 == '3' && p_n <= n_obs
        % w(x_{i,1},...,x_{i,p_n}) = 1/\sqrt{1+(x_i-c)^T \inv(\Sigma) (x_i-c)},
        % based on robust Mahalanobix distance

        if     options.choice_of_weight_function_in_robust_BD_estimation == 31
            % MCD estimator
            [hat_center, hat_Cov] = robust_center_scatter(X, 1);

        elseif options.choice_of_weight_function_in_robust_BD_estimation == 32
            % S-estimator
            [hat_center, hat_Cov] = robust_center_scatter(X, 2);
        end

        w_robust = zeros(n_obs, 1);
        for i = 1:n_obs
            X_i = X(i, :)'; % p_n*1 vector

            w_robust(i) = 1/sqrt(1 + (X_i - hat_center)' * pinv(hat_Cov) ...
                * (X_i - hat_center));
        end
    end
end
