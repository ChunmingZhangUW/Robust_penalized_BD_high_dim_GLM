
function [p_1_vector, p_2_vector, results] = robust_p_1_p_2_BD_for_hat_V_n(...
    I_loss, family, link, theta, y, options, index_robust_c, ...
    choice_rho_function, c_tune_constant)

%--------------------------------------------------------------------------
% Name    : robust_p_1_p_2_BD_for_hat_V_n.m
% Function: find pp_1(y;\theta) and pp_2(y;\theta)
%           used in weighted quadratic approximation
%           to robust Bregman Divergence \rho_q(y, F^{-1}(\theta))
% Model   : F(m(x))=beta_1 x_1 + ...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss    : deviance loss, exponential loss, quadratic loss,
%           quasi-likelihood
% Link F  : the canonical link function
% Used in : estimate the covariance matrix of the minimizer of
%           1/n*\sum_{i=1}^n \rho_q(Y_i, m(X_i)) w(X_i)
% Called  : q_1_q_2_BD_for_hat_V_n.m, robust_rho_function.m,
%           deri_G_1_function.m,
%           Binomial_pmf.m, Binomial_cdf.m,
%           Poisson_pmf.m,  Poisson_cdf.m,
%           hat_sigma2_of_Gaussian_data_in_robust_BD.m,
%           hat_phi_of_overdispersed_Poisson_data_in_robust_BD.m
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

if index_robust_c == 0 || c_tune_constant == inf  % \psi(r) = r

    [p_1_vector, p_2_vector] = q_1_q_2_BD_for_hat_V_n(...
        I_loss, family, link, theta, y, options);
    return
end

%---------------------- for robust-BD ------------------------

% index_robust_c == 1 && c_tune_constant ~= inf  % \psi(r) \ne r

n_obs = size(theta, 1);

%------------------------- compute mu --------------------------

if     family == 0 && strcmpi(link, 'iden')
    % Gaussian responses, identity link

    mu = theta;

    residual = y - mu;

    if     options.choice_of_estimating_sigma2_inside_iteration_robust_BD == 0
        % estimate sigma^2 outside iterations
        hat_sigma2 = options.hat_sigma2_outside_iteration;

    elseif options.choice_of_estimating_sigma2_inside_iteration_robust_BD == 1
        % estimate sigma^2 inside each iteration
        hat_sigma2 = hat_sigma2_of_Gaussian_data_in_robust_BD(...
            I_loss, family, link, theta, y, options);
    end

    results.V_original = hat_sigma2;  % V(u) = \sigma^2

    deri_V = 0;

elseif family == 1 && strcmpi(link, 'logit')
    % Bernoulli responses, logit link

    mu = 1./(1 + exp(-theta));

    residual = y - mu;

    results = [];

    deri_V = 1 - 2*mu;

elseif family == 12 && strcmpi(link, 'logit')
    % Binomial responses, logit link

    y_Bin = y(:, 1);
    N_Bin = y(:, 2);

    p_Bin = 1./(1 + exp(-theta));

    mu = N_Bin.*p_Bin;

    residual = y_Bin - mu;

    results = [];

    deri_V = 1 - 2*p_Bin;  % V'(u)

elseif family == 21 && strcmpi(link, 'log')
    % Poisson responses, log link

    theta(theta < -options.BD_C) = -options.BD_C;
    theta(theta >  options.BD_C) =  options.BD_C;

    mu = exp(theta);

    residual = y - mu;

    hat_Poisson_phi = hat_phi_of_overdispersed_Poisson_data_in_robust_BD(...
        I_loss, family, link, theta, y, options, c_tune_constant);
    Poisson_phi = hat_Poisson_phi;

    results.Poisson_phi = Poisson_phi;

    deri_V = Poisson_phi*1;
end

if     options.general_formula_robust_p_1_p_2_BD_for_hat_V_n == 1 ...
        || (family == 1  && strcmpi(link, 'logit') && c_tune_constant < 1) ...
        || (family == 12 && strcmpi(link, 'logit') && strcmpi(choice_rho_function, 'Tukey_biweight')) ...
        || (family == 21 && strcmpi(link, 'log')   && strcmpi(choice_rho_function, 'Tukey_biweight'))
    % Bernoulli responses, logit link, c_tune_constant < 1
    % Binomial  responses, logit link, Tukey-psi function
    % Poisson   responses, log   link, Tukey-psi function
    %----------------- general formulas for pp_1, pp_2 -------------------

    if     family == 0  % Gaussian responses
        V_original = hat_sigma2;  % V(u) = \sigma^2

        if strcmpi(link, 'iden')
            deri_1_F = 1;
            deri_2_F = 0;
        end

        if (I_loss == 1 || I_loss == 3) % deviance, quadratic
            deri_2_q = -2;
            deri_3_q = 0;
        end

    elseif family == 1  % Bernoulli responses
        one_minus_mu = 1 - mu;
        V_original = mu.*one_minus_mu;
        V_modified = max(options.zero_thres, V_original);

        if strcmpi(link, 'logit')
            deri_1_F = 1./V_modified;
            deri_2_F = -deri_V./(V_modified.^2);
        end

        if     I_loss == 1  % deviance
            deri_2_q = -2./V_modified;
            deri_3_q =  2*deri_V./(V_modified.^2);

        elseif I_loss == 2  % exponential
            deri_2_q = -1/2       ./(V_modified.^(3/2));
            deri_3_q =  3/4*deri_V./(V_modified.^(5/2));
        end

    elseif family == 12 % Binomial responses
        V_original = mu.*(1 - p_Bin);    % mu.*(1-p);
        V_modified = max(options.zero_thres, V_original);

        if strcmpi(link, 'logit')
            deri_1_F = 1./V_modified;
            deri_2_F = -deri_V./(V_modified.^2);
        end

        if I_loss == 1 % deviance loss
            deri_2_q = -2*deri_1_F;
            deri_3_q = -2*deri_2_F;
        end

    elseif family == 21 % Poisson responses
        V_original = Poisson_phi*mu;  % V(u) = \phi u
        mu_modified = max(options.zero_thres, mu);

        if strcmpi(link, 'log')
            deri_1_F =  1./mu_modified;
            deri_2_F = -1./(mu_modified.^2);
        end

        if I_loss == 211  % V(x)=phi x, (negative) quasi-likelihood
            deri_2_q = -deri_1_F;
            deri_3_q = -deri_2_F;
        end
    end

    sqrt_V_original = sqrt(V_original);
    if family == 0     % Gaussian responses
        sqrt_V_modified = sqrt_V_original;
    else
        sqrt_V_modified = max(options.zero_thres, sqrt_V_original);
    end

    r_y_mu = residual./sqrt_V_modified;  % Pearson_residual, r(y, mu)

    psi_r_y_mu       = robust_rho_function(r_y_mu, index_robust_c, ...
        choice_rho_function, c_tune_constant, 1);
    % \psi(r(y, mu))
    psi_prime_r_y_mu = robust_rho_function(r_y_mu, index_robust_c, ...
        choice_rho_function, c_tune_constant, 2);
    % \psi'(r(y, mu))

    deri_1_G_1 = deri_G_1_function(family, mu, y, results, options, ...
        index_robust_c, choice_rho_function, c_tune_constant, 1);
    % G_1'(mu)

    psi_c_minus_deri_1_G_1 = psi_r_y_mu - deri_1_G_1;
    % \psi(r(y, mu)) - G_1'(mu)
    A_0_mu = deri_2_q.*sqrt_V_original./deri_1_F; % A_0(mu)

    p_1_vector = (psi_c_minus_deri_1_G_1).*A_0_mu; % p_1(y; \theta)

    %---------------------------

    deri_1_psi_r_y_mu = - psi_prime_r_y_mu.*( 1 + ...
        r_y_mu.*deri_V./(2*sqrt_V_modified) )./sqrt_V_modified;
    % d\psi(r(y, mu))/d\mu
    deri_2_G_1 = deri_G_1_function(family, mu, y, results, options, ...
        index_robust_c, choice_rho_function, c_tune_constant, 2);
    % G_1''(mu)

    A_0_y_mu = (deri_1_psi_r_y_mu - deri_2_G_1).*A_0_mu./deri_1_F;
    % A_0(y, mu)
    A_1_mu = ...
        ( (deri_3_q.*sqrt_V_original ...
        + 1/2*deri_2_q.*deri_V./sqrt_V_modified).*deri_1_F ...
        - deri_2_q.*sqrt_V_original.*deri_2_F )./(deri_1_F.^3);
    % A_1(mu)

    p_2_vector = A_0_y_mu + (psi_c_minus_deri_1_G_1).*A_1_mu; % p_2(y; \theta)

    return

elseif options.general_formula_robust_p_1_p_2_BD_for_hat_V_n == 0
    %----------------- explicit formulas for pp_1, pp_2 -------------------

    if     family == 0 && strcmpi(link, 'iden')
        % Gaussian responses, identity link

        if (I_loss == 1 || I_loss == 3) % deviance, quadratic
            V_original = hat_sigma2;  % V(u) = \sigma^2
            sqrt_V_original = sqrt(V_original);

            if     strcmpi(choice_rho_function, 'Huber') % quadratic
                c_sqrt_V_original = c_tune_constant*sqrt_V_original;

                mu_1 = y - c_sqrt_V_original;
                mu_2 = y + c_sqrt_V_original;

                p_1_vector = ...
                    + 2*c_sqrt_V_original.*(mu >= mu_2) ...
                    - 2*residual         .*(mu_1 < mu & mu < mu_2) ...
                    - 2*c_sqrt_V_original.*(mu <= mu_1);

                %---------------------------
                if     options.Gaussian_robust_p_2_vector_in_robust_p_1_p_2_BD_for_hat_V_n == 0
                    % original p_2
                    p_2_vector = zeros(n_obs, 1);
                    p_2_vector(mu_1 < mu & mu < mu_2) = 2;

                elseif options.Gaussian_robust_p_2_vector_in_robust_p_1_p_2_BD_for_hat_V_n == 1
                    % p_2 is approximated by 2
                    p_2_vector = 2;

                elseif options.Gaussian_robust_p_2_vector_in_robust_p_1_p_2_BD_for_hat_V_n == 2
                    pearson_residual = residual/sqrt(hat_sigma2);       % r
                    psi_c_r = robust_rho_function(pearson_residual, ...
                        index_robust_c, choice_rho_function, c_tune_constant, 1);
                    % \psi(r)

                    p_2_vector = 2*psi_c_r./pearson_residual;
                    % 2*\psi'(r) is approximated by 2*\psi(r)/r
                end

            elseif strcmpi(choice_rho_function, 'Tukey_biweight') % quadratic
                r_y_mu = residual/sqrt_V_original;

                psi_c_r   = robust_rho_function(r_y_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 1);
                d_psi_c_r = robust_rho_function(r_y_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 2);

                p_1_vector = psi_c_r*(-2*sqrt_V_original);
                p_2_vector = 2*d_psi_c_r;
            end
        end

    elseif family == 1 && strcmpi(link, 'logit')
        % Bernoulli responses, logit link

        one_minus_mu = 1 - mu;
        V_original = mu.*one_minus_mu;
        sqrt_V_original = sqrt(V_original);

        c_square = c_tune_constant^2;
        C_1 = 1/(1 + c_square); C_2 = 1 - C_1;

        half_c = c_tune_constant/2;

        if     c_tune_constant < 1
            disp('!!!robust_p_1_p_2_BD_for_hat_V_n.m: N.A. return')
            return

        elseif c_tune_constant >= 1
            if     I_loss == 1    % deviance loss
                if     strcmpi(choice_rho_function, 'Huber') % deviance loss
                    c_sqrt_V_original = c_tune_constant*sqrt_V_original;

                    p_1_vector = ...
                        + (mu           + c_sqrt_V_original).*(  0 <= mu & mu <= C_1) ...
                        + (C_1 < mu & mu < C_2) ...
                        + (one_minus_mu + c_sqrt_V_original).*(C_2 <= mu & mu <= 1);

                    p_1_vector = -2*residual.*p_1_vector;

                    %---------------------------
                    I_1 = residual*half_c.*deri_V;
                    I_2 = -c_tune_constant*V_original;
                    I_12 = I_1 + I_2;

                    p_2_vector = ...
                        + ( +(    y - 2*mu).*sqrt_V_original + I_12 ) ...
                        .*(0 <= mu & mu < C_1) ...
                        - sqrt_V_original ...
                        .*(C_1 < mu & mu < C_2) ...
                        + ( -(1 + y - 2*mu).*sqrt_V_original + I_12 ) ...
                        .*(C_2 < mu & mu <= 1);

                    p_2_vector = -2*sqrt_V_original.*p_2_vector;

                elseif strcmpi(choice_rho_function, 'Tukey_biweight') % deviance loss
                    p_1_vector = zeros(n_obs, 1);
                    p_1_vector(0 <= mu & mu <= C_1) = ...
                        + mu          (0 <= mu & mu <= C_1) ...
                        .*( 1 - (1/c_square).*mu          (0 <= mu & mu <= C_1) ...
                        ./one_minus_mu(0 <= mu & mu <= C_1) ).^2;
                    p_1_vector(C_1 < mu & mu < C_2) = ...
                        + mu          (C_1 < mu & mu < C_2) ...
                        .*( 1 - (1/c_square).*mu          (C_1 < mu & mu < C_2) ...
                        ./one_minus_mu(C_1 < mu & mu < C_2) ).^2 ...
                        + one_minus_mu(C_1 < mu & mu < C_2) ...
                        .*( 1 - (1/c_square).*one_minus_mu(C_1 < mu & mu < C_2) ...
                        ./mu          (C_1 < mu & mu < C_2) ).^2;
                    p_1_vector(C_2 <= mu & mu <= 1) = ...
                        + one_minus_mu(C_2 <= mu & mu <= 1) ...
                        .*( 1 - (1/c_square).*one_minus_mu(C_2 <= mu & mu <= 1) ...
                        ./mu          (C_2 <= mu & mu <= 1) ).^2;

                    p_1_vector = -2*residual.*p_1_vector;

                    %---------------------------
                    p_2_vector = zeros(n_obs, 1);
                    p_2_vector(0 <= mu & mu <= C_1) = ...
                        + ( y(0 <= mu & mu <= C_1) - 2*mu(0 <= mu & mu <= C_1) ) ...
                        .*( 1 - (1/c_square)*mu(0 <= mu & mu <= C_1) ...
                        ./one_minus_mu(0 <= mu & mu <= C_1) ).^2 ...
                        + residual(0 <= mu & mu <= C_1).*mu(0 <= mu & mu <= C_1) ...
                        *2.*( 1 - (1/c_square)*mu(0 <= mu & mu <= C_1) ...
                        ./one_minus_mu(0 <= mu & mu <= C_1) ) ...
                        *(-1/c_square)./(one_minus_mu(0 <= mu & mu <= C_1).^2);
                    p_2_vector(C_1 < mu & mu < C_2) = ...
                        + ( y(C_1 < mu & mu < C_2) - 2*mu(C_1 < mu & mu < C_2) ) ...
                        .*( 1 - (1/c_square)*mu(C_1 < mu & mu < C_2) ...
                        ./one_minus_mu(C_1 < mu & mu < C_2) ).^2 ...
                        + residual(C_1 < mu & mu < C_2).*mu(C_1 < mu & mu < C_2) ...
                        *2.*( 1 - (1/c_square)*mu(C_1 < mu & mu < C_2) ...
                        ./one_minus_mu(C_1 < mu & mu < C_2) ) ...
                        *(-1/c_square)./(one_minus_mu(C_1 < mu & mu < C_2).^2) ...
                        - ( 1 + y(C_1 < mu & mu < C_2) - 2*mu(C_1 < mu & mu < C_2) ) ...
                        .*( 1 - (1/c_square)*one_minus_mu(C_1 < mu & mu < C_2) ...
                        ./mu(C_1 < mu & mu < C_2) ).^2 ...
                        + residual(C_1 < mu & mu < C_2).*one_minus_mu(C_1 < mu & mu < C_2) ...
                        *2.*( 1 - (1/c_square)*one_minus_mu(C_1 < mu & mu < C_2) ...
                        ./mu(C_1 < mu & mu < C_2) ) ...
                        *(-1/c_square)*(-1)./(mu(C_1 < mu & mu < C_2).^2);
                    p_2_vector(C_2 <= mu & mu <= 1) = ...
                        - ( 1 + y(C_2 <= mu & mu <= 1) - 2*mu(C_2 <= mu & mu <= 1) ) ...
                        .*( 1 - (1/c_square)*one_minus_mu(C_2 <= mu & mu <= 1) ...
                        ./mu(C_2 <= mu & mu <= 1) ).^2 ...
                        + residual(C_2 <= mu & mu <= 1).*one_minus_mu(C_2 <= mu & mu <= 1) ...
                        *2.*( 1 - (1/c_square)*one_minus_mu(C_2 <= mu & mu <= 1) ...
                        ./mu(C_2 <= mu & mu <= 1) ) ...
                        *(-1/c_square)*(-1)./(mu(C_2 <= mu & mu <= 1).^2);

                    p_2_vector = -2*V_original.*p_2_vector;
                end

            elseif I_loss == 2   % exponential loss
                if     strcmpi(choice_rho_function, 'Huber') % exponential loss
                    c_sqrt_V_original = c_tune_constant*sqrt_V_original;

                    p_1_vector = zeros(n_obs, 1);
                    p_1_vector(0 <= mu & mu <= C_1) = ...
                        sqrt( mu(0 <= mu & mu <= C_1)./one_minus_mu(0 <= mu & mu <= C_1) ) ...
                        + c_tune_constant;
                    p_1_vector(C_1 < mu & mu < C_2) = ...
                        1./sqrt_V_original(C_1 < mu & mu < C_2);
                    p_1_vector(C_2 <= mu & mu <= 1) = ...
                        sqrt( one_minus_mu(C_2 <= mu & mu <= 1)./mu(C_2 <= mu & mu <= 1) ) ...
                        + c_tune_constant;
                    p_1_vector = - residual/2.*p_1_vector;

                    %---------------------------
                    p_2_vector = zeros(n_obs, 1);
                    p_2_vector(0 <= mu & mu < C_1) = ...
                        (1/2)*sqrt_V_original(0 <= mu & mu < C_1) ...
                        .*( mu(0 <= mu & mu < C_1) ...
                        + c_sqrt_V_original(0 <= mu & mu < C_1) ) ...
                        - residual(0 <= mu & mu < C_1)/4 ...
                        .*sqrt(mu(0 <= mu & mu < C_1)./one_minus_mu(0 <= mu & mu < C_1));
                    p_2_vector(C_1 < mu & mu < C_2) = ...
                        (1/2)*sqrt_V_original(C_1 < mu & mu < C_2) ...
                        + residual(C_1 < mu & mu < C_2)/4 ...
                        .*deri_V(C_1 < mu & mu < C_2) ...
                        ./sqrt_V_original(C_1 < mu & mu < C_2);
                    p_2_vector(C_2 < mu & mu <= 1) = ...
                        (1/2)*sqrt_V_original(C_2 < mu & mu <= 1) ...
                        .*( one_minus_mu(C_2 < mu & mu <= 1) ...
                        + c_sqrt_V_original(C_2 < mu & mu <= 1) ) ...
                        + residual(C_2 < mu & mu <= 1)/4 ...
                        .*sqrt(one_minus_mu(C_2 < mu & mu <= 1)./mu(C_2 < mu & mu <= 1));

                elseif strcmpi(choice_rho_function, 'Tukey_biweight') % exponential loss
                    p_1_vector = zeros(n_obs, 1);
                    p_1_vector(0 <= mu & mu <= C_1) = ...
                        + 1./one_minus_mu(0 <= mu & mu <= C_1) ...
                        .*( 1 - (1/c_square).*mu          (0 <= mu & mu <= C_1) ...
                        ./one_minus_mu(0 <= mu & mu <= C_1) ).^2;
                    p_1_vector(C_1 < mu & mu < C_2) = ...
                        + 1./one_minus_mu(C_1 < mu & mu < C_2) ...
                        .*( 1 - (1/c_square).*mu          (C_1 < mu & mu < C_2) ...
                        ./one_minus_mu(C_1 < mu & mu < C_2) ).^2 ... %
                        + 1./mu          (C_1 < mu & mu < C_2) ...
                        .*( 1 - (1/c_square).*one_minus_mu(C_1 < mu & mu < C_2) ...
                        ./mu          (C_1 < mu & mu < C_2) ).^2;
                    p_1_vector(C_2 <= mu & mu <= 1) = ...
                        + 1./mu          (C_2 <= mu & mu <= 1) ...
                        .*( 1 - (1/c_square).*one_minus_mu(C_2 <= mu & mu <= 1) ...
                        ./mu          (C_2 <= mu & mu <= 1) ).^2;

                    p_1_vector = -residual/2.*sqrt_V_original.*p_1_vector;

                    %---------------------------
                    p_2_vector = zeros(n_obs, 1);

                    p_2_vector(0 <= mu & mu <= C_1) = ...
                        + ( (y(0 <= mu & mu <= C_1) - 2*mu(0 <= mu & mu <= C_1)) ...
                        .*( 1 - (1/c_square)*mu(0 <= mu & mu <= C_1) ...
                        ./one_minus_mu(0 <= mu & mu <= C_1) ).^2 ...
                        + residual(0 <= mu & mu <= C_1).*mu(0 <= mu & mu <= C_1)*2 ...
                        .*(1 - (1/c_square)*mu(0 <= mu & mu <= C_1) ...
                        ./one_minus_mu(0 <= mu & mu <= C_1)) ...
                        *(-1/c_square)./(one_minus_mu(0 <= mu & mu <= C_1).^2) ) ...
                        .*sqrt_V_original(0 <= mu & mu <= C_1) ...
                        - residual(0 <= mu & mu <= C_1) ...
                        .*( 1 - (1/c_square)*mu(0 <= mu & mu <= C_1) ...
                        ./one_minus_mu(0 <= mu & mu <= C_1) ).^2 ...
                        .*deri_V(0 <= mu & mu <= C_1)/2.*sqrt( mu(0 <= mu & mu <= C_1) ...
                        ./one_minus_mu(0 <= mu & mu <= C_1) );

                    p_2_vector(C_1 < mu & mu < C_2) = ...
                        + ( (y(C_1 < mu & mu < C_2) - 2*mu(C_1 < mu & mu < C_2)) ...
                        .*( 1 - (1/c_square)*mu(C_1 < mu & mu < C_2) ...
                        ./one_minus_mu(C_1 < mu & mu < C_2) ).^2 ...
                        + residual(C_1 < mu & mu < C_2).*mu(C_1 < mu & mu < C_2)*2 ...
                        .*(1 - (1/c_square)*mu(C_1 < mu & mu < C_2) ...
                        ./one_minus_mu(C_1 < mu & mu < C_2)) ...
                        *(-1/c_square)./(one_minus_mu(C_1 < mu & mu < C_2).^2) ) ...
                        .*sqrt_V_original(C_1 < mu & mu < C_2) ...
                        - residual(C_1 < mu & mu < C_2) ...
                        .*( 1 - (1/c_square)*mu(C_1 < mu & mu < C_2) ...
                        ./one_minus_mu(C_1 < mu & mu < C_2) ).^2 ...
                        .*deri_V(C_1 < mu & mu < C_2)/2.*sqrt( mu(C_1 < mu & mu < C_2) ...
                        ./one_minus_mu(C_1 < mu & mu < C_2) ) ...  %
                        + ( - ( 1 + y(C_1 < mu & mu < C_2) - 2*mu(C_1 < mu & mu < C_2)) ...
                        .*( 1 - (1/c_square)*one_minus_mu(C_1 < mu & mu < C_2) ...
                        ./mu(C_1 < mu & mu < C_2) ).^2 ...
                        + residual(C_1 < mu & mu < C_2) ...
                        .*one_minus_mu(C_1 < mu & mu < C_2)*2 ...
                        .*(1 - (1/c_square)*one_minus_mu(C_1 < mu & mu < C_2) ...
                        ./mu(C_1 < mu & mu < C_2)) ...
                        *(-1/c_square)*(-1)./(mu(C_1 < mu & mu < C_2).^2) ) ...
                        .*sqrt_V_original(C_1 < mu & mu < C_2) ...
                        - residual(C_1 < mu & mu < C_2) ...
                        .*( 1 - (1/c_square)*one_minus_mu(C_1 < mu & mu < C_2) ...
                        ./mu(C_1 < mu & mu < C_2) ).^2 ...
                        .*deri_V(C_1 < mu & mu < C_2)/2 ...
                        .*sqrt( one_minus_mu(C_1 < mu & mu < C_2) ...
                        ./mu(C_1 < mu & mu < C_2) );

                    p_2_vector(C_2 <= mu & mu <= 1) = ...
                        + ( - ( 1 + y(C_2 <= mu & mu <= 1) - 2*mu(C_2 <= mu & mu <= 1)) ...
                        .*( 1 - (1/c_square)*one_minus_mu(C_2 <= mu & mu <= 1) ...
                        ./mu(C_2 <= mu & mu <= 1) ).^2 ...
                        + residual(C_2 <= mu & mu <= 1) ...
                        .*one_minus_mu(C_2 <= mu & mu <= 1)*2 ...
                        .*(1 - (1/c_square)*one_minus_mu(C_2 <= mu & mu <= 1) ...
                        ./mu(C_2 <= mu & mu <= 1)) ...
                        *(-1/c_square)*(-1)./(mu(C_2 <= mu & mu <= 1).^2) ) ...
                        .*sqrt_V_original(C_2 <= mu & mu <= 1) ...
                        - residual(C_2 <= mu & mu <= 1) ...
                        .*( 1 - (1/c_square)*one_minus_mu(C_2 <= mu & mu <= 1) ...
                        ./mu(C_2 <= mu & mu <= 1) ).^2 ...
                        .*deri_V(C_2 <= mu & mu <= 1)/2 ...
                        .*sqrt( one_minus_mu(C_2 <= mu & mu <= 1) ...
                        ./mu(C_2 <= mu & mu <= 1) );

                    p_2_vector = -1/2*p_2_vector;
                end
            end
        end

    elseif family == 12 && strcmpi(link, 'logit')
        % Binomial responses, logit link

        if I_loss == 1 % deviance loss
            if     strcmpi(choice_rho_function, 'Huber') % deviance loss
                V_original = mu.*(1 - p_Bin);    % mu.*(1-p);
                sqrt_V_original = sqrt(V_original);

                c_sqrt_V_original = c_tune_constant*sqrt_V_original;

                J_1 = floor(mu - c_sqrt_V_original);
                J_2 =  ceil(mu + c_sqrt_V_original) - 1;

                N_minus_y = N_Bin - y_Bin;  % N - y
                c_square = c_tune_constant^2;

                mu_1 = ( (2*y_Bin + c_square).*N_Bin - c_tune_constant ...
                    *sqrt( N_Bin.*( 4*y_Bin.*N_minus_y + c_square*N_Bin ) ) ) ...
                    ./( 2*(N_Bin + c_square) );
                mu_2 = ( (2*y_Bin + c_square).*N_Bin + c_tune_constant ...
                    *sqrt( N_Bin.*( 4*y_Bin.*N_minus_y + c_square*N_Bin ) ) ) ...
                    ./( 2*(N_Bin + c_square) );

                %--------------------------------------------------------

                F_N_J_1                     = Binomial_cdf(J_1,   N_Bin,   p_Bin);
                F_N_minus_one_J_1_minus_one = Binomial_cdf(J_1-1, N_Bin-1, p_Bin);
                F_N_minus_two_J_1_minus_two = Binomial_cdf(J_1-2, N_Bin-2, p_Bin);

                F_N_J_2                     = Binomial_cdf(J_2,   N_Bin,   p_Bin);
                F_N_minus_one_J_2_minus_one = Binomial_cdf(J_2-1, N_Bin-1, p_Bin);
                F_N_minus_two_J_2_minus_two = Binomial_cdf(J_2-2, N_Bin-2, p_Bin);

                %--------------------

                I_11      = -c_tune_constant*F_N_J_1;

                I_11_star =  c_square       *F_N_J_1;

                I_13      = c_tune_constant*(1 - F_N_J_2);

                I_13_star = c_square       *(1 - F_N_J_2);

                I_12_0_0  = ...
                    + (F_N_minus_one_J_2_minus_one - F_N_minus_one_J_1_minus_one) ...
                    - (F_N_J_2                     - F_N_J_1);
                I_12_0    = mu.*I_12_0_0;

                %--------------------

                T_1_mu = F_N_J_2 - F_N_J_1;

                %--------------------

                I_21_times_V = mu*(-c_tune_constant).* ...
                    ( F_N_minus_one_J_1_minus_one - F_N_J_1 );

                I_23_times_V = mu*(c_tune_constant).* ...
                    ( F_N_J_2 - F_N_minus_one_J_2_minus_one );

                I_22_0 = ...
                    + mu.* (...
                    + (F_N_minus_two_J_2_minus_two - F_N_minus_two_J_1_minus_two) ...
                    - (F_N_minus_one_J_2_minus_one - F_N_minus_one_J_1_minus_one) ) ...
                    - p_Bin.* ...
                    (F_N_minus_two_J_2_minus_two - F_N_minus_two_J_1_minus_two) ...
                    + (F_N_minus_one_J_2_minus_one - F_N_minus_one_J_1_minus_one) ...
                    - I_12_0;
                %- mu.* ( ...
                %+ (F_N_minus_one_J_2_minus_one - F_N_minus_one_J_1_minus_one) ...
                %- (F_N_J_2 - F_N_J_1) );
                I_22_0 = max(I_22_0, 0); % I_22_0 >= 0

                %--------------------------------------------------------

                P_11_vector = ...
                    + 2*c_sqrt_V_original.*(mu >= mu_2) ...
                    - 2*(y_Bin - mu)     .*(mu_1 < mu & mu < mu_2) ...
                    - 2*c_sqrt_V_original.*(mu <= mu_1);

                P_12_vector = 2*( sqrt_V_original.*(I_11 + I_13) + I_12_0 );
                % 2\sqrt{V(u)} G_1'(u)

                p_1_vector  = P_11_vector + P_12_vector;

                %--------------------------------------------------------

                P_21_vector   = ...
                    + c_sqrt_V_original.*deri_V.*(mu > mu_2) ...
                    + 2*V_original             .*(mu_1 < mu & mu < mu_2) ...
                    - c_sqrt_V_original.*deri_V.*(mu < mu_1);

                P_22_1_vector = - V_original.*T_1_mu - deri_V/2.*I_12_0;
                P_22_1_vector = 2*P_22_1_vector;
                % 2\sqrt{V(u)} Delta_1(u) V(u)

                P_22_2_vector = sqrt_V_original.*(I_21_times_V + I_23_times_V) ...
                    + mu.*I_22_0;
                P_22_2_vector = max(P_22_2_vector, 0); % P_22_2_vector >= 0
                P_22_2_vector = 2*P_22_2_vector;
                % 2\sqrt{V(u)} Delta_2(u) V(u)

                P_22_vector   = P_22_1_vector + P_22_2_vector;
                % 2\sqrt{V(u)} G_1''(u) V(u)

                P_23_vector   = (P_12_vector/2).*deri_V;
                % \sqrt{V(u)} G_1'(u) V'(u)

                p_2_vector    = P_21_vector + P_22_vector + P_23_vector;

                %--------------------------------------------------------

                E_psi_c_square_times_V = (I_11_star + I_13_star).*V_original ...
                    + mu.*I_22_0;
                %E\{\psi^2(r(Y, u) | u\} V(u)
                results.E_psi_c_square_times_V = E_psi_c_square_times_V;

                deri_1_psi_c_r = (J_1 + 1 <= y_Bin & y_Bin <= J_2); % \psi'(r(y, u))
                results.deri_1_psi_c_r = deri_1_psi_c_r;

                results.P_12_vector   = P_12_vector;   % 2\sqrt{V(u)} G_1'(u)

                results.P_22_2_vector = P_22_2_vector; % 2\sqrt{V(u)} Delta_2(u) V(u)

            elseif strcmpi(choice_rho_function, 'Tukey_biweight') % deviance loss
                disp('!!!robust_p_1_p_2_BD_for_hat_V_n.m: N.A. return')
                return
            end
        end

    elseif family == 21 && strcmpi(link, 'log')
        % Poisson responses, log link

        if I_loss == 211  % V(x)=phi x, (negative) quasi-likelihood
            sqrt_Poisson_phi = sqrt(Poisson_phi);
            V_original = Poisson_phi*mu;  % V(u) = \phi u
            sqrt_V_original = sqrt(V_original);

            c_sqrt_V_original = c_tune_constant*sqrt_V_original;

            J_1 = floor(mu - c_sqrt_V_original);
            J_2 =  ceil(mu + c_sqrt_V_original) - 1;

            c_square = c_tune_constant^2;

            mu_1 = ( (2*y + Poisson_phi*c_square) ...
                - sqrt_Poisson_phi*c_tune_constant ...
                *sqrt(4*y + Poisson_phi*c_square) )/2;
            mu_2 = ( (2*y + Poisson_phi*c_square) ...
                + sqrt_Poisson_phi*c_tune_constant ...
                *sqrt(4*y + Poisson_phi*c_square) )/2;

            %--------------------------------------------------------

            if     Poisson_phi == 1 % genuine Poisson
                if     strcmpi(choice_rho_function, 'Huber') % quasi-likelihood
                    F_J_1 = Poisson_cdf(J_1, mu);
                    P_J_1 = Poisson_pmf(J_1, mu);

                    F_J_2 = Poisson_cdf(J_2, mu);
                    P_J_2 = Poisson_pmf(J_2, mu);

                    %--------------------

                    I_11      = -c_tune_constant*F_J_1;

                    I_11_star =  c_square       *F_J_1;

                    I_13      = c_tune_constant*(1 - F_J_2);

                    I_13_star = c_square       *(1 - F_J_2);

                    I_12_0_0  = P_J_1 - P_J_2;
                    I_12_0    = mu.*I_12_0_0;

                    %--------------------

                    T_1_mu = F_J_2 - F_J_1;  % T_1(u)

                    %--------------------

                    I_21 = c_tune_constant*P_J_1;

                    I_23 = c_tune_constant*P_J_2;

                    P_J_1_minus_1 = Poisson_pmf(J_1 - 1, mu);
                    P_J_2_minus_1 = Poisson_pmf(J_2 - 1, mu);
                    F_J_1_minus_1 = F_J_1 - P_J_1;
                    F_J_2_minus_1 = F_J_2 - P_J_2; %Poisson_cdf(J_2-1, mu);
                    I_22_0 =  ...
                        + mu.*( P_J_1_minus_1 - P_J_2_minus_1 ) ...
                        + ( F_J_2_minus_1 - F_J_1_minus_1 ) ...
                        - I_12_0;
                    I_22_0 = max(I_22_0, 0); % I_22_0 >= 0

                    %---------------------------

                    E_psi_c_square_times_V = (I_11_star + I_13_star).*V_original ...
                        + mu.*I_22_0;
                    %E\{\psi^2(r(Y, u) | u\} V(u)
                    results.E_psi_c_square_times_V = E_psi_c_square_times_V;

                elseif strcmpi(choice_rho_function, 'Tukey_biweight') % quasi-likelihood
                    disp('!!!robust_p_1_p_2_BD_for_hat_V_n.m: N.A. return')
                    return
                end

            elseif Poisson_phi ~= 1 % overdispersed Poisson
                if     strcmpi(choice_rho_function, 'Huber') % quasi-likelihood
                    delta = 1; % used in the derivative calculation

                    d_1 = zeros(n_obs, 1); d_2 = zeros(n_obs, 1);
                    for i = 1:n_obs
                        mu_i = mu(i);

                        dd_delta = abs( mu - (mu_i+delta) );
                        dd_delta = sort(dd_delta);
                        d_1(i) = dd_delta(2);

                        dd = abs( mu - mu_i );
                        dd = sort(dd);
                        d_2(i) = dd(2);
                    end
                    band_1 = 5*max(d_1);
                    band_2 = 5*max(d_2);

                    I_11   = zeros(n_obs, 1);
                    I_13   = zeros(n_obs, 1);
                    I_12_0 = zeros(n_obs, 1);

                    T_1_mu = zeros(n_obs, 1);

                    hat_F_J_1 = zeros(n_obs, 1);
                    hat_I_12_0      = zeros(n_obs, 1);
                    hat_F_J_2       = zeros(n_obs, 1);

                    hat_I_21   = zeros(n_obs, 1);
                    hat_I_23   = zeros(n_obs, 1);
                    hat_I_22_1 = zeros(n_obs, 1);
                    hat_I_22_2 = zeros(n_obs, 1);

                    %--------------------
                    for i = 1:n_obs
                        if band_1 <= eps
                            K_h_delta = ones(n_obs, 1);
                        else
                            argu_delta = ( mu - (mu_i+delta) )/band_1;
                            K_h_delta  = (1 - argu_delta.^2).*(abs(argu_delta) <= 1);
                        end
                        mean_K_h_delta = mean(K_h_delta);

                        if band_2 <= eps
                            K_h = ones(n_obs, 1);
                        else
                            argu = ( mu - mu_i )/band_2;
                            K_h  = (1 - argu.^2).*(abs(argu) <= 1);
                        end
                        mean_K_h = mean(K_h);

                        if     options.method_estimation == 1  % empirical
                            I_11(i)   = -c_tune_constant * mean( y <= J_1(i) ) ...
                                *(J_1(i) >= 0);
                            I_13(i)   =  c_tune_constant * mean( y >= J_2(i)+1 );
                            I_12_0(i) =  ...
                                mean( (y - mu_i) .* (J_1(i)+1 <= y & y <= J_2(i)) );

                            T_1_mu(i) = mean( J_1(i)+1 <= y & y <= J_2(i) );

                        elseif options.method_estimation == 2  % kernel
                            hat_F_J_1(i) = mean( K_h.*(y <= J_1(i)) )/mean_K_h;

                            hat_I_12_0(i) = ...
                                - mean( K_h .*( 0 <= y & y <= J_1(i) ).*(y - mu_i) ...
                                )/mean_K_h*(J_1(i) >= 0) ...
                                - mean( K_h .*( y >= J_2(i)+1 ).*(y - mu_i) ...
                                )/mean_K_h;

                            hat_F_J_2(i)       = mean( K_h.*( y <= J_2(i)     ) )/mean_K_h;
                        end

                        hat_I_21(i) = 1/delta *( ...
                            + mean( K_h_delta.*( y <= J_1(i) ) )...
                            /mean_K_h_delta ...
                            - mean( K_h      .*( y <= J_1(i) ) )...
                            /mean_K_h );

                        hat_I_23(i) = 1/delta *( ...
                            + mean( K_h_delta.*( y >= J_2(i)+1 ) )...
                            /mean_K_h_delta ...
                            - mean( K_h      .*( y >= J_2(i)+1 ) )...
                            /mean_K_h );

                        hat_I_22_1(i) = 1/delta *( ...
                            + mean( K_h_delta.*( y <= J_1(i) ).*(y - mu_i) )...
                            /mean_K_h_delta ...
                            - mean( K_h      .*( y <= J_1(i) ).*(y - mu_i) )...
                            /mean_K_h );

                        hat_I_22_2(i) = 1/delta *( ...
                            + mean( K_h_delta.*( y >= J_2(i)+1 ).*(y - mu_i) )...
                            /mean_K_h_delta ...
                            - mean( K_h      .*( y >= J_2(i)+1 ).*(y - mu_i) )...
                            /mean_K_h );
                    end
                    %--------------------

                    if options.method_estimation == 2  % kernel
                        I_11   = -c_tune_constant*hat_F_J_1;

                        I_13   =  c_tune_constant*(1 - hat_F_J_2);

                        I_12_0 =  hat_I_12_0;

                        T_1_mu = hat_F_J_2 - hat_F_J_1;
                    end

                    %--------------------

                    I_21   = -c_tune_constant    .*hat_I_21 ;

                    I_23   =  c_tune_constant    .*hat_I_23;

                    I_22_0 =  1 - hat_I_22_1 - hat_I_22_2;

                elseif strcmpi(choice_rho_function, 'Tukey_biweight') % quasi-likelihood
                    disp('!!!robust_p_1_p_2_BD_for_hat_V_n.m: N.A. return')
                    return
                end
            end % end of two cases of count data

            %--------------------------------------------------------

            P_11_vector = ...
                + c_sqrt_V_original.*(mu >= mu_2) ...
                - residual         .*(mu_1 < mu & mu < mu_2) ...
                - c_sqrt_V_original.*(mu <= mu_1);

            P_12_vector = sqrt_V_original.*(I_11 + I_13) + I_12_0;
            % \sqrt{V(u)} G_1'(u)

            p_1_vector  = P_11_vector + P_12_vector;

            %--------------------------------------------------------

            P_21_vector   = ...
                + c_sqrt_V_original/2.*(mu > mu_2) ...
                + mu                 .*(mu_1 < mu & mu < mu_2) ...
                - c_sqrt_V_original/2.*(mu < mu_1);

            P_22_1_vector = - mu.*T_1_mu - I_12_0/2;
            %\sqrt{V(u)} Delta_1(u) u

            P_22_2_vector = sqrt_V_original.*(I_21 + I_23).*mu + mu.*I_22_0;
            P_22_2_vector = max(P_22_2_vector, 0); % P_22_2_vector >= 0
            %\sqrt{V(u)} Delta_2(u) u

            P_22_vector   = P_22_1_vector + P_22_2_vector;

            P_23_vector   = P_12_vector/2;
            % \sqrt{V(u)} G_1'(u)/2

            p_2_vector    = P_21_vector + P_22_vector + P_23_vector;

            %--------------------------------------------------------

            deri_1_psi_c_r = (mu_1 < mu & mu < mu_2); % \psi'(r(y, u))
            results.deri_1_psi_c_r = deri_1_psi_c_r;

            results.P_12_vector   = P_12_vector;   % \sqrt{V(u)} G_1'(u)

            results.P_22_2_vector = P_22_2_vector; % \sqrt{V(u)} Delta_2(u) u
        end

        % else
        %     sqrt_V_original = sqrt(V_original);
        %     sqrt_V_modified = max(options.zero_thres, sqrt_V_original);
        %
        %     r_y_mu = residual./sqrt_V_modified;  % Pearson_residual
        %
        %     psi      = robust_rho_function(r_y_mu, index_robust_c, ...
        %         choice_rho_function, c_tune_constant, 1);
        %     % 1st derivative of Huber func.
        %     psi_prime_r_y_mu = robust_rho_function(r_y_mu, index_robust_c, ...
        %         choice_rho_function, c_tune_constant, 2);
        %     % 2nd derivative of Huber func.
        %
        %     deri_1_G_1 = G_1(mu, sqrt_V_original, sqrt_V_modified, family, options, ...
        %         index_robust_c, c_tune_constant, 1);
        %     psi_c_minus_deri_1_G_1 = psi - deri_1_G_1;
        %
        %     deri_2_G_1 = G_1(mu, sqrt_V_original, sqrt_V_modified, family, options, ...
        %         index_robust_c, c_tune_constant, 2);
        %
        %     p_1_vector = (psi_c_minus_deri_1_G_1).*deri_2_q.*sqrt_V_original./deri_1_F;
        %
        %     %---------------------------
        %     deri_r = -(V + residual.*deri_V/2)./(sqrt_V_modified.*V_modified);
        %     deri_psi_c_multiply_deri_r = psi_prime_r_y_mu.*deri_r; % \psi'(r)r'(u)
        %     f_1 = (deri_psi_c_multiply_deri_r - deri_2_G_1).*deri_2_q.*sqrt_V_original ...
        %         ./deri_1_F;
        %     f_2 = (psi_c_minus_deri_1_G_1).*...
        %         ( (deri_3_q.*sqrt_V_original+deri_2_q.*deri_V./sqrt_V_modified/2) ...
        %         .*deri_1_F - deri_2_q.*sqrt_V_original.*deri_2_F )./(deri_1_F.^2);
        %
        %     p_2_vector = (f_1 + f_2)./deri_1_F;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p_2_vector = p_2_vector + (p_2_vector == 0)*options.delta + ...
%     sign(p_2_vector)*options.delta;
% This makes problems to
% robust_LARS_parametric_BD_estimate_SCAD_LLA_general.m.

%p_2_vector = max(p_2_vector, options.delta);  % set to be >= options.delta

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if any(isempty(p_1_vector)) == 1
    disp(' !!!robust_p_1_p_2_BD_for_hat_V_n.m: some estimate of p_1_vector = [ ]!!!');
end

if any(isnan(p_1_vector))   == 1
    disp(' !!!robust_p_1_p_2_BD_for_hat_V_n.m: some estimate of p_1_vector = NaN!!!');
end

if any(isinf(p_1_vector))   == 1
    disp(' !!!robust_p_1_p_2_BD_for_hat_V_n.m: some estimate of p_1_vector = Inf!!!');
end
%-------------

if any(isempty(p_2_vector)) == 1
    disp(' !!!robust_p_1_p_2_BD_for_hat_V_n.m: some estimate of p_2_vector = [ ]!!!');
end

if any(isnan(p_2_vector))   == 1
    disp(' !!!robust_p_1_p_2_BD_for_hat_V_n.m: some estimate of p_2_vector = NaN!!!');
end

if any(isinf(p_2_vector))   == 1
    disp(' !!!robust_p_1_p_2_BD_for_hat_V_n.m: some estimate of p_2_vector = Inf!!!');
end