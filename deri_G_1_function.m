
function f = deri_G_1_function(family, mu, y, results, options, ...
    index_robust_c, choice_rho_function, c_tune_constant, deri)

%--------------------------------------------------------------------------
% Name    : deri_G_1_function.m
% Function: compute G_1'(mu) and G_1''(mu)
%--------------------------------------------------------------------------
% <Input>
% index_robust_c:  0: \psi(r)  =  r;
%                  1: \psi(r) \ne r;
%choice_rho_function: 1 Huber \rho function; 2 Tukey biweight function
% c_tune_constant: constant used in \psi(r)
%--------------------------------------------------------------------------
% <Output>
%--------------------------------------------------------------------------

if deri == 0 || deri >= 3
    disp(' !!!deri_G_1_function.m: derivative is not equal to 1 or 2, return')
    return
end

n_obs = size(mu, 1);
f = zeros(n_obs, 1);  % initialize

%--------------------- for classical-BD ----------------------

if index_robust_c == 0 || c_tune_constant == inf  % \psi(r) = r

    f = zeros(n_obs, 1);
    return
end

%---------------------- for robust-BD ------------------------

% index_robust_c == 1 && c_tune_constant ~= inf  % \psi(r) \ne r

if     family == 0    % Gaussian responses
    if     deri == 1
        f = zeros(n_obs, 1);

    elseif deri == 2
        f = zeros(n_obs, 1);
    end

elseif family == 1    % Bernoulli responses
    one_minus_mu = 1 - mu;
    V_original = mu.*one_minus_mu;
    sqrt_V_original = sqrt(V_original);

    c_square = c_tune_constant^2;
    C_1 = 1/(1 + c_square); C_2 = 1 - C_1;

    if     c_tune_constant < 1
        sqrt_V_modified = max(options.zero_thres, sqrt_V_original);

        if     deri == 1
            f = 0;
            for j = 0:1
                P_j_mu = (mu.^j).*(one_minus_mu.^(1-j)); % P(Y=j)=[1-mu, mu]

                r_j_mu = (j - mu)./sqrt_V_modified;      % r(j, mu)

                psi_r_j_mu = robust_rho_function(r_j_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 1);  % \psi(r(j, mu))

                f = f + psi_r_j_mu.*P_j_mu;
            end
            % G_1'(mu)=\sum_{j=0}^1 \psi(r(j, mu))*P(Y=j)

        elseif deri == 2
            deri_V = 1 - 2*mu; % V'(mu)

            f = 0;
            for j = 0:1
                P_j_mu = (mu.^j).*(one_minus_mu.^(1-j)); % P(Y=j)=[1-mu, mu]

                r_j_mu = (j - mu)./sqrt_V_modified;      % r(j, mu)

                psi_r_j_mu = robust_rho_function(r_j_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 1);  % \psi(r(j, mu))

                psi_prime_r_j_mu = robust_rho_function(r_j_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 2);  % \psi'(r(j, mu))

                deri_1_psi_r_j_mu = - psi_prime_r_j_mu.*( 1 + ...
                    r_j_mu.*deri_V./(2*sqrt_V_modified) )./sqrt_V_modified;
                % d\psi(r(j, mu))/d\mu

                deri_1_P_j_mu = - (j==0) + (j==1);
                % d P_j_mu/d\mu

                f = f ...
                    + deri_1_psi_r_j_mu.*P_j_mu ...
                    + psi_r_j_mu       .*deri_1_P_j_mu;
            end
            % G_1''(mu) = \sum_{j=0}^1 [d\psi(r(j, mu))/d\mu]*P(Y=j) ...
            %            +\sum_{j=0}^1 \psi(r(j, mu))*[d P(Y=j)/d\mu]
        end

    elseif c_tune_constant >= 1
        if     strcmpi(choice_rho_function, 'Huber')
            if     deri == 1
                f = ...
                    +(-sqrt_V_original + c_tune_constant*mu) ...
                    .*(0   <= mu & mu <= C_1) ...
                    +(+sqrt_V_original - c_tune_constant*one_minus_mu) ...
                    .*(C_2 <= mu & mu <= 1);

            elseif deri == 2
                deri_V = 1 - 2*mu;
                sqrt_V_modified = max(options.zero_thres, sqrt_V_original);

                f = ...
                    +(-deri_V./(2*sqrt_V_modified) + c_tune_constant) ...
                    .*(0   <= mu & mu <  C_1) ...
                    +( deri_V./(2*sqrt_V_modified) + c_tune_constant) ...
                    .*(C_2 <  mu & mu <= 1);
            end

        elseif strcmpi(choice_rho_function, 'Tukey_biweight')
            if     deri == 1
                f(0 <= mu & mu <= C_1) = ...
                    -sqrt_V_original(0 <= mu & mu <= C_1) ...
                    .*( 1 - (1/c_square)*mu(0 <= mu & mu <= C_1) ...
                    ./(one_minus_mu(0 <= mu & mu <= C_1)) ).^2;
                f(C_1 < mu & mu < C_2) = ...
                    -sqrt_V_original(C_1 < mu & mu < C_2) ...
                    .*( 1 - (1/c_square)*mu(C_1 < mu & mu < C_2) ...
                    ./(one_minus_mu(C_1 < mu & mu < C_2)) ).^2 ...
                    +sqrt_V_original(C_1 < mu & mu < C_2) ...
                    .*( 1 - (1/c_square)*one_minus_mu(C_1 < mu & mu < C_2) ...
                    ./(mu(C_1 < mu & mu < C_2)) ).^2;
                f(C_2 <= mu & mu <= 1) = ...
                    +sqrt_V_original(C_2 <= mu & mu <= 1) ...
                    .*( 1 - (1/c_square)*one_minus_mu(C_2 <= mu & mu <= 1) ...
                    ./(mu(C_2 <= mu & mu <= 1)) ).^2;

            elseif deri == 2
                sqrt_V_modified = max(options.zero_thres, sqrt_V_original);
                deri_V = 1 - 2*mu;

                f(0 <= mu & mu <= C_1) = ...
                    - ( deri_V(0 <= mu & mu <= C_1) ...
                    ./( 2*sqrt_V_modified(0 <= mu & mu <= C_1) ) ...
                    .*( 1 - (1/c_square)*mu(0 <= mu & mu <= C_1) ...
                    ./(one_minus_mu(0 <= mu & mu <= C_1)) ).^2 ...
                    + sqrt_V_original(0 <= mu & mu <= C_1) ...
                    *2.*( 1 - (1/c_square)*mu(0 <= mu & mu <= C_1) ...
                    ./(one_minus_mu(0 <= mu & mu <= C_1)) ) ...
                    *(-1/c_square)./(one_minus_mu(0 <= mu & mu <= C_1).^2) );
                f(C_1 < mu & mu < C_2) = ...
                    - ( deri_V(C_1 < mu & mu < C_2) ...
                    ./(2*sqrt_V_modified(C_1 < mu & mu < C_2)) ...
                    .*( 1 - (1/c_square)*mu(C_1 < mu & mu < C_2) ...
                    ./(one_minus_mu(C_1 < mu & mu < C_2)) ).^2 ...
                    + sqrt_V_original(C_1 < mu & mu < C_2) ...
                    *2.*( 1 - (1/c_square)*mu(C_1 < mu & mu < C_2) ...
                    ./(one_minus_mu(C_1 < mu & mu < C_2)) ) ...
                    *(-1/c_square)./(one_minus_mu(C_1 < mu & mu < C_2).^2) ) ...
                    + ( deri_V(C_1 < mu & mu < C_2) ...
                    ./(2*sqrt_V_modified(C_1 < mu & mu < C_2)) ...
                    .*( 1 - ( 1/c_square)*one_minus_mu(C_1 < mu & mu < C_2) ...
                    ./(mu(C_1 < mu & mu < C_2)) ).^2 ...
                    + sqrt_V_original(C_1 < mu & mu < C_2) ...
                    *2.*( 1 - (1/c_square)*one_minus_mu(C_1 < mu & mu < C_2) ...
                    ./(mu(C_1 < mu & mu < C_2)) ) ...
                    *(1/c_square)./(mu(C_1 < mu & mu < C_2).^2) );
                f(C_2 <= mu & mu <= 1) = ...
                    + ( deri_V(C_2 <= mu & mu <= 1) ...
                    ./(2*sqrt_V_modified(C_2 <= mu & mu <= 1)) ...
                    .*( 1 - (1/c_square)*one_minus_mu(C_2 <= mu & mu <= 1) ...
                    ./(mu(C_2 <= mu & mu <= 1)) ).^2 ...
                    + sqrt_V_original(C_2 <= mu & mu <= 1) ...
                    *2.*( 1 - (1/c_square)*one_minus_mu(C_2 <= mu & mu <= 1) ...
                    ./(mu(C_2 <= mu & mu <= 1)) ) ...
                    *(1/c_square)./(mu(C_2 <= mu & mu <= 1).^2) );
            end
        end
    end

elseif family == 12    % Binomial responses
    N_Bin = y(:, 2);

    p_Bin = mu./N_Bin;

    V_original = mu.*(1 - p_Bin);    % mu.*(1-p);
    sqrt_V_original = sqrt(V_original);
    sqrt_V_modified = max(options.zero_thres, sqrt_V_original);

    c_sqrt_V_original = c_tune_constant*sqrt_V_original;

    J_1 = floor(mu - c_sqrt_V_original);
    J_2 =  ceil(mu + c_sqrt_V_original) - 1;

    F_N_J_1                     = Binomial_cdf(J_1,   N_Bin,   p_Bin);
    F_N_minus_one_J_1_minus_one = Binomial_cdf(J_1-1, N_Bin-1, p_Bin);
    F_N_minus_two_J_1_minus_two = Binomial_cdf(J_1-2, N_Bin-2, p_Bin);

    F_N_J_2                     = Binomial_cdf(J_2,   N_Bin,   p_Bin);
    F_N_minus_one_J_2_minus_one = Binomial_cdf(J_2-1, N_Bin-1, p_Bin);
    F_N_minus_two_J_2_minus_two = Binomial_cdf(J_2-2, N_Bin-2, p_Bin);

    if     strcmpi(choice_rho_function, 'Huber')
        if     deri == 1
            %--------------------

            I_11      = -c_tune_constant*F_N_J_1;

            I_13      =  c_tune_constant*(1 - F_N_J_2);

            I_12_0_0  = ...
                + (F_N_minus_one_J_2_minus_one - F_N_minus_one_J_1_minus_one) ...
                - (F_N_J_2                     - F_N_J_1);
            I_12_0    = mu.*I_12_0_0;
            I_12      = I_12_0./sqrt_V_modified;

            %--------------------

            f = I_11 + I_13 + I_12;

        elseif deri == 2
            deri_V = 1 - 2*p_Bin;  % V'(mu)
            V_modified = max(options.zero_thres, V_original);

            %--------------------------------

            T_1_mu = F_N_J_2 - F_N_J_1;

            I_12_0_0  = ...
                + (F_N_J_1 - F_N_minus_one_J_1_minus_one) ...
                - (F_N_J_2       - F_N_minus_one_J_2_minus_one);
            I_12_0    = mu.*I_12_0_0;
            T_2_mu = I_12_0;

            Delta_1_mu = ( - T_1_mu - deri_V./(2*V_modified).*T_2_mu ) ...
                ./sqrt_V_modified;

            %--------------------------------

            I_21_times_V = mu*(-c_tune_constant).* ...
                ( F_N_minus_one_J_1_minus_one - F_N_J_1 );
            I_21 = I_21_times_V./V_modified;

            I_23_times_V = mu*(c_tune_constant).* ...
                ( F_N_J_2 - F_N_minus_one_J_2_minus_one );
            I_23 = I_23_times_V./V_modified;

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
            I_22 = (mu.*I_22_0./V_modified)./sqrt_V_modified;

            Delta_2_mu = I_21 + I_23 + I_22;

            %--------------------------------

            f = Delta_1_mu + Delta_2_mu;
        end

    elseif strcmpi(choice_rho_function, 'Tukey_biweight')
        if     deri == 1
            f = 0;
            for j = 0:N_Bin
                P_j_mu = Binomial_pmf(j, N_Bin, p_Bin);   % P(Y=j)

                r_j_mu = (j - mu)./sqrt_V_modified;       % r(j, mu)

                psi_r_j_mu = robust_rho_function(r_j_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 1);  % \psi(r(j, mu))

                f = f + psi_r_j_mu.*P_j_mu;
            end
            % G_1'(mu)=\sum_{j=0}^N \psi(r(j, mu))*P(Y=j)

        elseif deri == 2
            deri_V = 1 - 2*p_Bin;  % V'(mu)
            V_modified = max(options.zero_thres, V_original);

            f = 0;
            for j = 0:N_Bin
                P_j_mu = Binomial_pmf(j, N_Bin, p_Bin);   % P(Y=j)

                r_j_mu = (j - mu)./sqrt_V_modified;       % r(j, mu)

                psi_r_j_mu = robust_rho_function(r_j_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 1);  % \psi(r(j, mu))

                psi_prime_r_j_mu = robust_rho_function(r_j_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 2);  % \psi'(r(j, mu))

                deri_1_psi_r_j_mu = - psi_prime_r_j_mu.*( 1 + ...
                    r_j_mu.*deri_V./(2*sqrt_V_modified) )./sqrt_V_modified;
                % d\psi(r(j, mu))/d\mu

                deri_1_P_j_mu = P_j_mu.*(j - mu)./V_modified;
                % d P_j_mu/d\mu

                f = f ...
                    + deri_1_psi_r_j_mu.*P_j_mu ...
                    + psi_r_j_mu       .*deri_1_P_j_mu;
            end
            % G_1''(mu) = \sum_{j=0}^N [d\psi(r(j, mu))/d\mu]*P(Y=j) ...
            %            +\sum_{j=0}^N \psi(r(j, mu))*[d P(Y=j)/d\mu]
        end
    end

elseif family == 21    % Poisson_quasi
    Poisson_phi = results.Poisson_phi;

    V_original = Poisson_phi*mu;  % V(mu) = \phi mu
    sqrt_V_original = sqrt(V_original);
    sqrt_V_modified = max(options.zero_thres, sqrt_V_original);

    c_sqrt_V_original = c_tune_constant*sqrt_V_original;

    if     strcmpi(choice_rho_function, 'Huber')
        J_1 = floor(mu - c_sqrt_V_original);
        J_2 =  ceil(mu + c_sqrt_V_original) - 1;

        F_J_1 = Poisson_cdf(J_1, mu);
        P_J_1 = Poisson_pmf(J_1, mu);

        F_J_2 = Poisson_cdf(J_2, mu);
        P_J_2 = Poisson_pmf(J_2, mu);

        if     deri == 1
            %--------------------

            I_11      = -c_tune_constant*F_J_1;

            I_13      =  c_tune_constant*(1 - F_J_2);

            I_12_0_0  = P_J_1 - P_J_2;
            I_12_0    = mu.*I_12_0_0;
            I_12 = I_12_0./sqrt_V_modified;

            %--------------------

            f = I_11 + I_13 + I_12;

        elseif deri == 2
            V_modified = max(options.zero_thres, V_original);

            %--------------------------------
            T_1_mu = F_J_2 - F_J_1;  % T_1(mu)

            I_12_0_0  = P_J_1 - P_J_2;
            I_12_0    = mu.*I_12_0_0;
            T_2_mu = I_12_0;

            Delta_1_mu = ( - T_1_mu - Poisson_phi./(2*V_modified).*T_2_mu ) ...
                ./sqrt_V_modified;

            %--------------------------------

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
            I_22 = I_22_0./sqrt_V_modified;

            Delta_2_mu = I_21 + I_23 + I_22;

            %--------------------------------

            f = Delta_1_mu + Delta_2_mu;
        end

    elseif strcmpi(choice_rho_function, 'Tukey_biweight')
        N_Poi = floor(mu + c_sqrt_V_original); % upper bound in the summations

        if     deri == 1
            f = 0;
            for j = 0:N_Poi
                P_j_mu = Poisson_pmf(j, mu);         % P(Y=j)

                r_j_mu = (j - mu)./sqrt_V_modified;  % r(j, mu)

                psi_r_j_mu = robust_rho_function(r_j_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 1);  % \psi(r(j, mu))

                f = f + psi_r_j_mu.*P_j_mu;
            end
            % G_1'(mu)=\sum_{j=0}^N \psi(r(j, mu))*P(Y=j)

        elseif deri == 2
            deri_V = Poisson_phi;  % V'(mu)
            V_modified = max(options.zero_thres, V_original);

            f = 0;
            for j = 0:N_Poi
                P_j_mu = Poisson_pmf(j, mu);         % P(Y=j)

                r_j_mu = (j - mu)./sqrt_V_modified;  % r(j, mu)

                psi_r_j_mu = robust_rho_function(r_j_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 1);  % \psi(r(j, mu))

                psi_prime_r_j_mu = robust_rho_function(r_j_mu, index_robust_c, ...
                    choice_rho_function, c_tune_constant, 2);  % \psi'(r(j, mu))

                deri_1_psi_r_j_mu = - psi_prime_r_j_mu.*( 1 + ...
                    r_j_mu.*deri_V./(2*sqrt_V_modified) )./sqrt_V_modified;
                % d\psi(r(j, mu))/d\mu

                deri_1_P_j_mu = P_j_mu.*(j - mu)./V_modified;
                % d P_j_mu/d\mu

                f = f ...
                    + deri_1_psi_r_j_mu.*P_j_mu ...
                    + psi_r_j_mu       .*deri_1_P_j_mu;
            end
            % G_1''(mu) = \sum_{j=0}^N [d\psi(r(j, mu))/d\mu]*P(Y=j) ...
            %            +\sum_{j=0}^N \psi(r(j, mu))*[d P(Y=j)/d\mu]
        end
    end
end

if any(isempty(f)) == 1
    disp([' !!!deri_G_1_function.m: some estimate of f = [ ]!!!', ...
        ', deri = ', num2str(deri)]);
end

if any(isnan(f)) == 1
    disp([' !!!deri_G_1_function.m: some estimate of f = NaN!!!', ...
        ', deri = ', num2str(deri)]);
end

if any(isinf(f)) == 1
    disp([' !!!deri_G_1_function.m: some estimate of f = Inf!!!', ...
        ', deri = ', num2str(deri)]);
end
