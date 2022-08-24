
function f = robust_true_pp(I_loss, family, link, theta, y, deri, ...
    options, index_robust_c, choice_rho_function, c_tune_constant)

%--------------------------------------------------------------------------
% Name    : robust_true_pp.m
% Function: compute pp_j(y;\theta), i.e. values and derivatives of
%           robust Bregman Divergence
%           \rho_q(y, F^{-1}(\theta)) w.r.t. \theta
% Loss    : deviance loss, exponential loss, quadratic loss,
%           quasi-likelihood
% Called  : true_qq.m, deri_1_G.m, robust_rho_function.m,
%           robust_p_1_p_2_BD_for_hat_V_n.m, p_star_function.m,
%           G_function.m, hat_sigma2_of_Gaussian_data_in_robust_BD.m
%--------------------------------------------------------------------------
% <Input>
% I_loss : choice of loss function:
%          1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%          211 (V(x) = phi x)
% family : 0 (Gaussian); 1 (Bernoulli), 12 (Binomial); 21 (Poisson_quasi)
%  link  : type of link function, 'iden', 'logit', 'log'
%  theta : scalar or vector
%    y   : scalar or vector with the same size of theta
%   deri : 0, 1, 2 derivatives of the loss
% options: set of parameters
% index_robust_c:  0: \psi(r)  =  r;
%                  1: \psi(r) \ne r;
%choice_rho_function: 1 Huber \rho function; 2 Tukey biweight function
% c_tune_constant: constant used in \psi(r)
%--------------------------------------------------------------------------
% <Output>
%    f   : values and derivatives of the robust-BD, same size of y
%--------------------------------------------------------------------------

%--------------------- for classical-BD ----------------------

if index_robust_c == 0 || c_tune_constant == inf  % \psi(r) = r

    f = true_qq(I_loss, family, link, theta, y, deri, options);
    return
end

%---------------------- for robust-BD ------------------------

% index_robust_c == 1 && c_tune_constant ~= inf  % \psi(r) \ne r

n_obs = size(theta, 1);

%------------------------- compute mu --------------------------

if     family == 0 && strcmpi(link, 'iden')
    % Gaussian responses, identity link

    mu = theta;
    %
    %     residual = y - mu;
    %
    %     %---------------------------------------------------------------------
    %     if     options.choice_of_Gaussian_sigma2_in_robust_p_1_p_2_BD_for_hat_V_n == 1
    %         % MAD (robust method)
    %
    %         median_residual= median(residual);
    %
    %         hat_sigma = median(abs(residual-median_residual))/0.6745;
    %         hat_sigma2_MAD = hat_sigma^2;
    %
    %         hat_sigma2 = max(hat_sigma2_MAD);
    %
    %     elseif options.choice_of_Gaussian_sigma2_in_robust_p_1_p_2_BD_for_hat_V_n == 2
    %         % max(MSE, MAD)
    %
    %         median_residual= median(residual); square_residual = residual.^2;
    %
    %         hat_sigma = median(abs(residual-median_residual))/0.6745;
    %         hat_sigma2_MAD = hat_sigma^2;
    %
    %         hat_sigma2_MSE = mean(square_residual);
    %
    %         hat_sigma2_Median_SE = median(square_residual);
    %
    %         hat_sigma2 = max([hat_sigma2_Median_SE, hat_sigma2_MSE, ...
    %             hat_sigma2_MAD]);
    %
    %     elseif options.choice_of_Gaussian_sigma2_in_robust_p_1_p_2_BD_for_hat_V_n == 3  % specified
    %
    %         disp(' --> enter robust_true_pp.m: for Gaussian with specified sigma2')
    %         hat_sigma2 = options.specified_robust_p_1_p_2_BD_Gaussian_sigma2;
    %
    %     else
    %         disp('!!! robust_true_pp.m: hat_Gaussian_sigma2 is undefined.')
    %         return
    %     end
    %     %---------------------------------------------------------------------

    if     options.choice_of_estimating_sigma2_inside_iteration_robust_BD == 0
        % estimate sigma^2 outside iterations
        hat_sigma2 = options.hat_sigma2_outside_iteration;

    elseif options.choice_of_estimating_sigma2_inside_iteration_robust_BD == 1
        % estimate sigma^2 inside each iteration
        hat_sigma2 = hat_sigma2_of_Gaussian_data_in_robust_BD(...
            I_loss, family, link, theta, y, options);
    end

    results.V_original = hat_sigma2;  % V(u) = \sigma^2

elseif family == 1 && strcmpi(link, 'logit')
    % Bernoulli responses, logit link

    mu = 1./(1 + exp(-theta));

    results = [];

elseif family == 12 && strcmpi(link, 'logit')
    % Binomial responses, logit link

    y_Bin = y(:, 1);
    N_Bin = y(:, 2);

    p_Bin = 1./(1 + exp(-theta));

    mu = N_Bin.*p_Bin;

    results = [];

elseif family == 21 && strcmpi(link, 'log')
    % Poisson responses, log link

    theta(theta < -options.BD_C) = -options.BD_C;
    theta(theta >  options.BD_C) =  options.BD_C;

    mu = exp(theta);
    %
    %     residual = y - mu;
    %
    %     %---------------------------------------------------------------------
    %     Poisson_a_UB = min(( (10^15-mu)/c_tune_constant ).^2./mu);
    %     square_residual = residual.^2;
    %
    %     if     options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n == 1
    %
    %         hat_Poisson_phi = min([...
    %             mean(square_residual.*mu)/mean(mu.^2), ...
    %             Poisson_a_UB, options.Poisson_variance_a_UB]);
    %
    %     elseif options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n == 2
    %
    %         hat_Poisson_phi = min([max([...
    %             median(square_residual./mu), ...
    %             mean(square_residual./mu) ...
    %             ]), ...
    %             Poisson_a_UB, options.Poisson_variance_a_UB]);
    %
    %     elseif options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n == 3  % specified
    %         % applicable to the genuine Poisson responses with phi = 1
    %
    %         %disp(' --> enter robust_true_pp.m: for Poisson with specified phi')
    %         hat_Poisson_phi = options.specified_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n;
    %
    %     else
    %         disp('!!! robust_true_pp.m: hat_Poisson_phi is undefined.')
    %         return
    %     end
    %
    %     Poisson_phi = hat_Poisson_phi;
    %     %---------------------------------------------------------------------

    hat_Poisson_phi = hat_phi_of_overdispersed_Poisson_data_in_robust_BD(...
        I_loss, family, link, theta, y, options, c_tune_constant);
    Poisson_phi = hat_Poisson_phi;

    results.Poisson_phi = Poisson_phi;
end

if     options.general_formula_robust_true_pp == 1
    %----------------- general formula for robust_true_pp -------------------

    if     deri == 0
        if family == 12    % Binomial responses
            LL_p_star = y_Bin;

        else
            LL_p_star = y; % lower limit in the intergration for p^*
        end
        LL_G = LL_p_star;  % lower limit in the intergration for G

        UL = mu;           % upper limit in the intergration

        p_star = zeros(n_obs, 1);
        G_mu   = zeros(n_obs, 1);
        for i = 1:n_obs
            p_star(i) = integral( ...
                @(x) deri_1_p_star(family, I_loss, x, index_robust_c, ...
                choice_rho_function, c_tune_constant, y(i, :), results, ...
                options), LL_p_star(i), UL(i));  % p^*(y, mu)

            G_mu(i) = integral( ...
                @(x) deri_1_G     (family, I_loss, x, index_robust_c, ...
                choice_rho_function, c_tune_constant, y(i, :), results, ...
                options), LL_G(i),      UL(i));  % G(mu)
        end

        p_0_vector = p_star - G_mu;

        f = p_0_vector;

    elseif deri == 1 || deri == 2

        [p_1_vector, p_2_vector] = robust_p_1_p_2_BD_for_hat_V_n(I_loss, ...
            family, link, theta, y, options, index_robust_c, ...
            choice_rho_function, c_tune_constant);

        if     deri == 1

            f = p_1_vector;

        elseif deri == 2

            f = p_2_vector;
        end
    end

elseif options.general_formula_robust_true_pp == 0
    %----------------- explicit formula for robust_true_pp -------------------

    if     deri == 0

        p_star_y_mu = p_star_function(I_loss, family, theta, y, mu, results, ...
            options, index_robust_c, choice_rho_function, c_tune_constant);

        if     family == 12  % Binomial responses

            p_star_y_y  = p_star_function(I_loss, family, theta, y, y_Bin, results, ...
                options, index_robust_c, choice_rho_function, c_tune_constant);

        elseif family ~= 12  % non-Binomial responses

            p_star_y_y  = p_star_function(I_loss, family, theta, y, y,     results, ...
                options, index_robust_c, choice_rho_function, c_tune_constant);
        end

        p_star = p_star_y_mu - p_star_y_y;

        %----------------------------------------------------------

        G_mu = G_function(I_loss, family, y, mu, results, ...
            options, index_robust_c, choice_rho_function, c_tune_constant);

        p_0_vector = p_star - G_mu;

        f = p_0_vector;

    elseif deri == 1 || deri == 2

        [p_1_vector, p_2_vector] = robust_p_1_p_2_BD_for_hat_V_n(I_loss, ...
            family, link, theta, y, options, index_robust_c, ...
            choice_rho_function, c_tune_constant);

        if     deri == 1

            f = p_1_vector;

        elseif deri == 2

            f = p_2_vector;
        end
    end

end


if any(isempty(f)) == 1
    disp([' !!!robust_true_pp.m: some estimate of f = [ ]!!!', ', I_loss=', ...
        num2str(I_loss), ', deri = ', num2str(deri)]);
end

if any(isnan(f)) == 1
    disp([' !!!robust_true_pp.m: some estimate of f = NaN!!!', ', I_loss=', ...
        num2str(I_loss), ', deri = ', num2str(deri)]);
end

if any(isinf(f)) == 1
    disp([' !!!robust_true_pp.m: some estimate of f = Inf!!!', ', I_loss=', ...
        num2str(I_loss), ', deri = ', num2str(deri)]);
end
