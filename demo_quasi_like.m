%
% Name    : demo_quasi_like.m
% Function: Penalized quasi-likelihood estimation for Poisson data
% Penalty : L1 and SCAD; weighted-L1
% Data    : simulated from overdispersed Poission responses
% Called  : generate_Poisson.m,
%           Poisson_initial_value_beta_0_in_robust_GLM_BD.m,
%           robust_GLM_BD_CD_penalized_parameter_estimate_SCAD.m,
%           robust_GLM_BD_CD_penalized_parameter_estimate.m;
%           GLM_BD_initial_weights_PCR.m,
%           robust_GLM_BD_parameter_estimate.m,
%           robust_GLM_BD_tune_lambda.m, robust_GLM_BD_tune_lambda_SCAD.m,
%           true_qq.m
%--------------------------------------------------------------------------
clear;
close all;

time = clock;
fprintf('%d-%d-%d %2d:%2d:%2.0f \n', time(1:6));

%------------------------------------------------------------
seeds_set = [1234 1234 1234 1234 1234 1234 1234 1234 1234];

n_common    = input('input sample size n (100; 200) = '); %200; %100;
disp(' ')
num_obs_set = [n_common, n_common]; % common n for different dimension p_n
p_n_set     = [ 50, 500]; %[ 50, 200];
num_sim     = input('input # (100; 500) of simulations = '); %100*5; %100*2;

%------------------------------------------------------------
norm_EE = 2; % |\hat \beta-\beta|_{norm_EE}

family = 21;  % Poisson responses
link = 'log';
I_loss = 211; % V(x)=phi x

rho = 0.2; %input(' input rho (0, 0.2) = ');
example = 1;

phi = 2; %input(' input phi (1: for genuine Poisson; 2: for overdispersed Poisson) = ');
options.true_Poisson_phi = phi; % true phi, used for generating data
% dispersion parameter
options.true_Poisson = (phi == 1);

misclass_Ber = [];
% misclass_Ber =  1 % compute misclassification rate for Bernoulli responses;
% misclass_Ber ~= 1 % compute other deviance
%------------------------------------------------------------

disp(' ')

options.CD_LARS = 1; % method (1: CD) for penalized estimation

lambda_method = 0; % input(' input method (0: tuning set; 1: CV) for selecting lambda = ');

disp(' ');

options.linear_log_grid = 2; % 1 for linear grid; 2 for log-scale grid
options.choice_initial_value = 1;
options.maxit = 40;
options.thresh_1 = 1e-3;     % for |b_1-b_0|
options.thresh_2 = 1e-3;     % for |b_1|
options.delta = eps; %1e-10;
options.eps = eps;           % for adding increase
options.max_abs_hat_beta = 1e4; % max(abs(hat_beta)) in parametric_BD_estimate.m
options.zero_thres = 1e-10;
options.SCAD_a = 3.7;
if     options.CD_LARS == 1
    options.CD_SCAD_LLA = 1;   % 1: CD   algorithm with SCAD penalty uses LLA
    options.CD_active_set = true;
end
if family == 21  % Poisson responses
    options.BD_C = 10;
    options.Poisson_variance_a_UB = 4;
end

study = input(' input study (1: without contamination; 2: with contamination) = ');
if     study == 1
    choice_contamination = 0;
    A = [];

elseif study == 2
    choice_contamination = 1;
    A = 50; %input(' input A (50; 100) = ');
end

index_robust = [0 0];  % initialize
index_robust(1) = input(' input index_robust (0: without; 1: with) robust for Y with \psi_c = ');
index_robust(2) = input(' input index_robust (0: without; 1: with) robust for X with weight = ');
disp(' ')

%------------------------------------ robust-BD ---------------
if     options.true_Poisson == 1  % genuine Poisson
    options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n = 3;
    options.specified_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n = 1;

elseif options.true_Poisson == 0  % overdispersed Poisson
    options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n = 3;
    %input([' input options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n \n',...
    %'  (1: Ratio_1; 2: Ratio_2; 3: specified) = ']);

    if     options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n ~= 3
        % estimated
        options.specified_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n = [];

    elseif options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n == 3
        % specified
        options.specified_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n = 1;
        %input(...
        %    ' input options.specified_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n = ');
    end
end
if options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n ~= 3 || ...
        options.specified_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n ~= 1
    options.method_estimation = input(...
        ' input options.method_estimation (1: empirical; 2: kernel) = ');
end
%-----------------------------------------------------------------
disp(' ')

number_choice_rho_function = 1;
%input(' input number_choice_rho_function (1: Huber; 2: Tukey) = ');
options.general_formula_robust_true_pp = 0;
%input('options.general_formula_robust_true_pp (1: general formula; 0: explicit formula) = ');
options.general_formula_robust_p_1_p_2_BD_for_hat_V_n = 0;
%input('options.general_formula_robust_p_1_p_2_BD_for_hat_V_n (1: general formula; 0: explicit formula) = ');
disp(' ');

specified_c_tune_constant_Huber = 1.345;
specified_c_tune_constant_Tukey = 4.685;

if     number_choice_rho_function == 1
    options.choice_rho_function = 'Huber';
    options.c_tune_constant     = specified_c_tune_constant_Huber;

elseif number_choice_rho_function == 2
    options.choice_rho_function = 'Tukey_biweight';
    options.c_tune_constant     = specified_c_tune_constant_Tukey;
end
if     options.c_tune_constant == 1.345
    c_tune_constant_character = '1345';
elseif options.c_tune_constant == 4.685
    c_tune_constant_character = '4685';
else
    c_tune_constant_character = input(' convert c_tune_constant (0.8) to character (''08'') = ');
    % used in formin the output figure file in .pdf format
end


if     index_robust(1) == 0  % \psci_c(r)  =  r
    choice_rho_function = [];
    c_tune_constant     = inf;

elseif index_robust(1) == 1  % \psci_c(r) \ne r
    choice_rho_function = options.choice_rho_function;
    c_tune_constant     = options.c_tune_constant;
end

options.choice_of_weight_function_in_robust_BD_estimation = 11;
%input([' options.choice_of_weight_function_in_robust_BD_estimation \n',...
%'  (10: s_n=1; 11: MAD est. of s_n; 12: Rousseeuw and Croux est. of s_n; 2: sqrt(1-hii)) = ']);
disp(' ')

options.GLM_BD_parameter_estimate = 21;
%input(' options.GLM_BD_parameter_estimate (direct: 1; CD: 21, 22, 23) = ');
options.robust_GLM_BD_parameter_estimate = 21;
%input(' options.robust_GLM_BD_parameter_estimate (direct: 1; CD: 21, 22, 23) = ');
if     options.CD_LARS == 1
    options.GLM_BD_CD_penalized_parameter_estimate = 21;
    %input(' input options.GLM_BD_CD_penalized_parameter_estimate (CD: 21, 22) = ');
    options.robust_GLM_BD_CD_penalized_parameter_estimate = 21;
    %input(' input options.robust_GLM_BD_CD_penalized_parameter_estimate (CD: 21, 22) = ');
end
disp(' ')

for cases = 1:length(num_obs_set)
    disp([' cases=', num2str(cases)]);
    num_obs = num_obs_set(cases);
    p_n = p_n_set(cases);     K = p_n + 1;

    sparse_beta = [5/2 2 2]'; % true non-zero coefficients
    true_beta = zeros(K, 1);
    true_beta(1:length(sparse_beta)) = sparse_beta;
    loc_zero = find(true_beta == 0);    % locations where true_beta = 0
    loc_nonzero = find(true_beta ~= 0); % locations where true_beta \ne 0
    s = sum(true_beta(2:end) ~= 0);     % s = 5

    penalty_set = (2 : length(true_beta));

    %######################## Part 1 ############################
    rng(seeds_set(cases), 'twister');

    n_grid_lambda = 11;
    %grid_A = log10(2^(-10)); grid_B = log10(2^(10));
    %grid_B = log10(2^(-10)); grid_A = log10(2^(10));
    grid_B = log10(2^(-10)); grid_A = log10(2^(0));

    kappa_grid = 2.^(-4 : 1 : 4); n_grid_kappa = length(kappa_grid);

    n_tune = 5000;
    % n_tune is the size of tuning sample for calculating test error

    j_to_display = length(sparse_beta) + 3;
    bias_beta_PCR = zeros(num_sim, j_to_display, 1);

    beta_error_SCAD   = zeros(num_sim, 1);
    beta_error_L1     = zeros(num_sim, 1);
    beta_error_CR     = zeros(num_sim, 1);
    beta_error_PCR    = zeros(num_sim, 1);
    beta_error_ORACLE = zeros(num_sim, 1);

    beta_MSE_SCAD   = zeros(num_sim, 1);   beta_MSE_L1  = zeros(num_sim, 1);
    beta_MSE_CR     = zeros(num_sim, 1);   beta_MSE_PCR = zeros(num_sim, 1);
    beta_MSE_ORACLE = zeros(num_sim, 1);

    ME_SCAD   = zeros(num_sim, 1);   ME_L1  = zeros(num_sim, 1);
    ME_CR     = zeros(num_sim, 1);   ME_PCR = zeros(num_sim, 1);
    ME_ORACLE = zeros(num_sim, 1);

    RME_SCAD   = zeros(num_sim, 1);   RME_L1 = zeros(num_sim, 1);
    RME_CR     = zeros(num_sim, 1);   RME_PCR = zeros(num_sim, 1);
    RME_ORACLE = zeros(num_sim, 1);

    TE_true = zeros(num_sim, 1);   TE_SCAD   = zeros(num_sim, 1);
    TE_L1   = zeros(num_sim, 1);   TE_CR     = zeros(num_sim, 1);
    TE_PCR  = zeros(num_sim, 1);   TE_ORACLE = zeros(num_sim, 1);

    number_CZ_SCAD   = zeros(num_sim, 1);  number_CNZ_SCAD   = zeros(num_sim, 1);
    number_CZ_L1     = zeros(num_sim, 1);  number_CNZ_L1     = zeros(num_sim, 1);
    number_CZ_CR     = zeros(num_sim, 1);  number_CNZ_CR     = zeros(num_sim, 1);
    number_CZ_PCR    = zeros(num_sim, 1);  number_CNZ_PCR    = zeros(num_sim, 1);
    number_CZ_ORACLE = zeros(num_sim, 1);  number_CNZ_ORACLE = zeros(num_sim, 1);
    %
    tic
    %
    for sim = 1:num_sim

        disp([' sim=', num2str(sim)]);

        %-------------- generate data set -------------------------

        converge = false;
        while converge == false
            [design_matrix, y_vector] = ...
                generate_Poisson(num_obs, p_n, true_beta, family, link, ...
                options.true_Poisson_phi, rho, study, example, A);

            [Beta_0, converge] = ...
                Poisson_initial_value_beta_0_in_robust_GLM_BD(...
                I_loss, family, link, ...
                design_matrix(:,loc_nonzero), y_vector, options);

            if converge == false
                disp('  within the while loop, converge = false; regenerate data')
            end
        end
        if converge == false
            disp('converge = false')
        end

        if lambda_method == 0 % tuning set
            [design_test_matrix, y_test_vector] = generate_Poisson(...
                num_obs, p_n, true_beta, family, link, ...
                options.true_Poisson_phi, rho, study, example, A);
            % 1 test set, num_obs*(p_n+1) for selecting lambda
        end

        [design_tune_matrix, y_tune_vector] = ...
            generate_Poisson(n_tune, p_n, true_beta, ...
            family, link, options.true_Poisson_phi, rho, study, example, A);
        % 1 tuning set, n_tune*(p_n+1) for computing ME & TE

        %------------------------------------------------------------

        w_robust = weight_function_in_robust_BD_estimation(...
            design_matrix, options, index_robust(2));

        I_loss = 211; % V(x)=x

        beta_0 = zeros(K, 1);
        beta_0(loc_nonzero) = Beta_0;

        % L1 penalty
        weight_pen = ones(length(penalty_set), 1); % penalty weights

        if     lambda_method == 0 % tuning set
            [lambda, Test_Error_L1, Hat_beta, I] = robust_GLM_BD_tune_lambda(...
                I_loss, family, link, ...
                design_matrix, y_vector, beta_0, ...
                design_test_matrix, y_test_vector, weight_pen, penalty_set, ...
                n_grid_lambda, grid_A, grid_B, misclass_Ber, options, ...
                index_robust(1), choice_rho_function, c_tune_constant, w_robust);
            hat_beta_L1 = Hat_beta(:, I);
        end % if lambda_method == 0 or 1

        % Oracle estimator
        beta_0_0 = beta_0;

        w_robust_ORACLE = weight_function_in_robust_BD_estimation(...
            design_matrix(:,loc_nonzero), options, index_robust(2));

        hat_beta_ORACLE = zeros(K, 1);
        [hat_beta_ORACLE(loc_nonzero), converge] = ...
            robust_GLM_BD_parameter_estimate( I_loss, family, link, ...
            design_matrix(:,loc_nonzero), y_vector, ...
            beta_0_0(loc_nonzero), options, index_robust(1), ...
            choice_rho_function, c_tune_constant, w_robust_ORACLE);

        if converge == false
            disp(' !!! fails to converge')
        end

        % SCAD penalty
        if     lambda_method == 0 % tuning set
            [lambda, Test_Error_SCAD, Hat_beta, I] = ...
                robust_GLM_BD_tune_lambda_SCAD( I_loss, family, link, ...
                design_matrix, y_vector, beta_0, ...
                design_test_matrix, y_test_vector, penalty_set, ...
                n_grid_lambda, grid_A, grid_B, misclass_Ber, options, ...
                index_robust(1), choice_rho_function, c_tune_constant, w_robust);
            hat_beta_SCAD =  Hat_beta(:, I);

        end % if lambda_method == 0 or 1

        % Weighted-L1 penalty, MR weight selection
        beta_init_CR = GLM_BD_initial_weights_PCR(I_loss, family, link, ...
            design_matrix(:, penalty_set), y_vector, 0, options);
        weight_pen = 1./abs(beta_init_CR); % penalty weights

        if     lambda_method == 0 % tuning set
            [lambda, Test_Error_CR, Hat_beta, I] = robust_GLM_BD_tune_lambda(...
                I_loss, family, link, ...
                design_matrix, y_vector, beta_0, ...
                design_test_matrix, y_test_vector, weight_pen, penalty_set, ...
                n_grid_lambda, grid_A, grid_B, misclass_Ber, options, ...
                index_robust(1), choice_rho_function, c_tune_constant, w_robust);
            hat_beta_CR = Hat_beta(:, I);
        end % if lambda_method == 0 or 1

        % Weighted-L1 penalty, PMR weight selection
        Hat3D_beta_PCR = zeros(K, n_grid_lambda, n_grid_kappa);
        V_row = zeros(1, n_grid_kappa); I_row = zeros(1, n_grid_kappa);
        Weight_Pen = zeros(length(penalty_set), n_grid_kappa);
        Lambda = zeros(n_grid_kappa, 1);

        for k = 1:n_grid_kappa
            kappa = kappa_grid(k);
            beta_init_PCR = GLM_BD_initial_weights_PCR(I_loss, family, link, ...
                design_matrix(:, penalty_set), y_vector, kappa, options);
            weight_pen = 1./abs(beta_init_PCR); % penalty weights

            if     lambda_method == 0 % tuning set
                [lambda, min_TE, Hat_beta, I] = robust_GLM_BD_tune_lambda(...
                    I_loss, family, link, design_matrix, y_vector, beta_0, ...
                    design_test_matrix, y_test_vector, ...
                    weight_pen, penalty_set, n_grid_lambda, grid_A, ...
                    grid_B, misclass_Ber, options, index_robust(1), ...
                    choice_rho_function, c_tune_constant, w_robust);

                I_row(k) = I;
                Hat3D_beta_PCR(:, :, k) = Hat_beta;
            end % if lambda_method == 0 or 1

            V_row(k) = min_TE;
        end % for k = 1:n_grid_kappa
        [CVE_PCR, kappa_opt] = min(V_row); % minimum CV error

        if     lambda_method == 0 % tuning set
            hat_beta_PCR = Hat3D_beta_PCR(:, I_row(kappa_opt), kappa_opt);
        end % if lambda_method == 0 or 1

        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        %--------------- regression estimation --------------------------

        bias_beta_PCR(sim, :, 1) = ...
            hat_beta_PCR(1:j_to_display) - true_beta(1:j_to_display);
        % Weighted-L1 PMR

        beta_error_SCAD  (sim, 1) = norm(hat_beta_SCAD   - true_beta, norm_EE);
        beta_error_L1    (sim, 1) = norm(hat_beta_L1     - true_beta, norm_EE);
        beta_error_CR    (sim, 1) = norm(hat_beta_CR     - true_beta, norm_EE);
        beta_error_PCR   (sim, 1) = norm(hat_beta_PCR    - true_beta, norm_EE);
        beta_error_ORACLE(sim, 1) = norm(hat_beta_ORACLE - true_beta, norm_EE);

        beta_MSE_SCAD  (sim, 1) = mean((hat_beta_SCAD   - true_beta).^2);
        beta_MSE_L1    (sim, 1) = mean((hat_beta_L1     - true_beta).^2);
        beta_MSE_CR    (sim, 1) = mean((hat_beta_CR     - true_beta).^2);
        beta_MSE_PCR   (sim, 1) = mean((hat_beta_PCR    - true_beta).^2);
        beta_MSE_ORACLE(sim, 1) = mean((hat_beta_ORACLE - true_beta).^2);

        %------------------ Test error --------------------

        true_theta       = design_tune_matrix*true_beta;
        hat_theta_SCAD   = design_tune_matrix*hat_beta_SCAD;
        hat_theta_L1     = design_tune_matrix*hat_beta_L1;
        hat_theta_CR     = design_tune_matrix*hat_beta_CR;
        hat_theta_PCR    = design_tune_matrix*hat_beta_PCR;
        hat_theta_ORACLE = design_tune_matrix*hat_beta_ORACLE;

        if     family == 0
            true_m       = true_theta;
            hat_m_SCAD   = hat_theta_SCAD;
            hat_m_L1     = hat_theta_L1;
            hat_m_CR     = hat_theta_CR;
            hat_m_PCR    = hat_theta_PCR;
            hat_m_ORACLE = hat_theta_ORACLE;

        elseif family == 1
            true_m       = 1./(1+exp(-true_theta));
            hat_m_SCAD   = 1./(1+exp(-hat_theta_SCAD));
            hat_m_L1     = 1./(1+exp(-hat_theta_L1));
            hat_m_CR     = 1./(1+exp(-hat_theta_CR));
            hat_m_PCR    = 1./(1+exp(-hat_theta_PCR));
            hat_m_ORACLE = 1./(1+exp(-hat_theta_ORACLE));

        elseif family == 21
            true_m       = exp(true_theta);
            hat_m_SCAD   = exp(hat_theta_SCAD);
            hat_m_L1     = exp(hat_theta_L1);
            hat_m_CR     = exp(hat_theta_CR);
            hat_m_PCR    = exp(hat_theta_PCR);
            hat_m_ORACLE = exp(hat_theta_ORACLE);
        end

        ME_SCAD  (sim, 1) = mean((hat_m_SCAD   - true_m).^2);
        ME_L1    (sim, 1) = mean((hat_m_L1     - true_m).^2);
        ME_CR    (sim, 1) = mean((hat_m_CR     - true_m).^2);
        ME_PCR   (sim, 1) = mean((hat_m_PCR    - true_m).^2);
        ME_ORACLE(sim, 1) = mean((hat_m_ORACLE - true_m).^2);

        RME_SCAD  (sim, 1) = ME_SCAD  (sim, 1)/ME_L1(sim, 1);
        RME_L1    (sim, 1) = ME_L1    (sim, 1)/ME_L1(sim, 1);
        RME_CR    (sim, 1) = ME_CR    (sim, 1)/ME_L1(sim, 1);
        RME_PCR   (sim, 1) = ME_PCR   (sim, 1)/ME_L1(sim, 1);
        RME_ORACLE(sim, 1) = ME_ORACLE(sim, 1)/ME_L1(sim, 1);

        predict_BD_true   = true_qq(I_loss, family, link, true_theta,       ...
            y_tune_vector, 0, options);
        predict_BD_SCAD   = true_qq(I_loss, family, link, hat_theta_SCAD,   ...
            y_tune_vector, 0, options);
        predict_BD_L1     = true_qq(I_loss, family, link, hat_theta_L1,     ...
            y_tune_vector, 0, options);
        predict_BD_CR     = true_qq(I_loss, family, link, hat_theta_CR,     ...
            y_tune_vector, 0, options);
        predict_BD_PCR    = true_qq(I_loss, family, link, hat_theta_PCR,    ...
            y_tune_vector, 0, options);
        predict_BD_ORACLE = true_qq(I_loss, family, link, hat_theta_ORACLE, ...
            y_tune_vector, 0, options);

        TE_true  (sim, 1) = mean(predict_BD_true);
        TE_SCAD  (sim, 1) = mean(predict_BD_SCAD);
        TE_L1    (sim, 1) = mean(predict_BD_L1);
        TE_CR    (sim, 1) = mean(predict_BD_CR);
        TE_PCR   (sim, 1) = mean(predict_BD_PCR);
        TE_ORACLE(sim, 1) = mean(predict_BD_ORACLE);

        %--------------- variable selection --------------------------

        number_CZ_SCAD  (sim, 1) = sum(hat_beta_SCAD  (loc_zero) == 0);
        number_CZ_L1    (sim, 1) = sum(hat_beta_L1    (loc_zero) == 0);
        number_CZ_CR    (sim, 1) = sum(hat_beta_CR    (loc_zero) == 0);
        number_CZ_PCR   (sim, 1) = sum(hat_beta_PCR   (loc_zero) == 0);
        number_CZ_ORACLE(sim, 1) = sum(hat_beta_ORACLE(loc_zero) == 0);

        number_CNZ_SCAD  (sim, 1) = sum(hat_beta_SCAD  (loc_nonzero) ~= 0);
        number_CNZ_L1    (sim, 1) = sum(hat_beta_L1    (loc_nonzero) ~= 0);
        number_CNZ_CR    (sim, 1) = sum(hat_beta_CR    (loc_nonzero) ~= 0);
        number_CNZ_PCR   (sim, 1) = sum(hat_beta_PCR   (loc_nonzero) ~= 0);
        number_CNZ_ORACLE(sim, 1) = sum(hat_beta_ORACLE(loc_nonzero) ~= 0);

    end % for sim = 1:num_sim
    t = toc;

    %%%%%%%%%%%%%%%%%%%% Output results %%%%%%%%%%%%%%%%%%%%%%%
    format short;
    disp([' phi = ', num2str(phi), ', rho = ', num2str(rho), ...
        ', example = ', num2str(example), ', sparse_beta = ', num2str(sparse_beta')])
    disp(['-- Case = ', num2str(cases), ...
        ', n_grid_lambda = ', num2str(n_grid_lambda), ...
        ', grid_A = ', num2str(grid_A), ', grid_B = ', num2str(grid_B)]);

    disp([' num_obs = ', num2str(num_obs), ', p_n = ', num2str(p_n), ...
        ', num_sim = ', num2str(num_sim), ...
        ', options.CD_LARS = ', num2str(options.CD_LARS), ...
        ', lambda_method = ', num2str(lambda_method), ...
        ', misclass_Ber = ', num2str(misclass_Ber), ...
        ', loss = ', num2str(I_loss), ', Poisson']);
    disp([' choice_contamination = ', num2str(choice_contamination), ...
        ', index_robust = [', num2str(index_robust), ...
        '], c_tune_constant = [', num2str(c_tune_constant), ']'])
    disp(['  TE_true = ', num2str(mean(TE_true))]);
    %disp('----------- regression estimates --------------');
    result_RME = [
        mean(RME_SCAD  (:, 1))
        mean(RME_L1    (:, 1))
        mean(RME_CR    (:, 1))
        mean(RME_PCR   (:, 1))
        mean(RME_ORACLE(:, 1))
        ];

    result_EE = [
        mean(beta_error_SCAD  (:, 1))  std(beta_error_SCAD  (:, 1));
        mean(beta_error_L1    (:, 1))  std(beta_error_L1    (:, 1));
        mean(beta_error_CR    (:, 1))  std(beta_error_CR    (:, 1));
        mean(beta_error_PCR   (:, 1))  std(beta_error_PCR   (:, 1));
        mean(beta_error_ORACLE(:, 1))  std(beta_error_ORACLE(:, 1))
        ];

    result_MSE = [
        mean(beta_MSE_SCAD  (:, 1))
        mean(beta_MSE_L1    (:, 1))
        mean(beta_MSE_CR    (:, 1))
        mean(beta_MSE_PCR   (:, 1))
        mean(beta_MSE_ORACLE(:, 1))
        ];

    result_ME = [
        mean(ME_SCAD  (:, 1))
        mean(ME_L1    (:, 1))
        mean(ME_CR    (:, 1))
        mean(ME_PCR   (:, 1))
        mean(ME_ORACLE(:, 1))
        ];

    %disp('----------- test errors --------------');
    result_TE = [
        mean(TE_SCAD  (:, 1))
        mean(TE_L1    (:, 1))
        mean(TE_CR    (:, 1))
        mean(TE_PCR   (:, 1))
        mean(TE_ORACLE(:, 1))
        ];

    %disp('------------ variable selection ---------------');
    result_VS = [
        mean(number_CZ_SCAD   (:, 1)) std(number_CZ_SCAD   (:, 1)) ...
        mean(number_CNZ_SCAD  (:, 1)) std(number_CNZ_SCAD  (:, 1));
        mean(number_CZ_L1     (:, 1)) std(number_CZ_L1     (:, 1)) ...
        mean(number_CNZ_L1    (:, 1)) std(number_CNZ_L1    (:, 1));
        mean(number_CZ_CR     (:, 1)) std(number_CZ_CR     (:, 1)) ...
        mean(number_CNZ_CR    (:, 1)) std(number_CNZ_CR    (:, 1));
        mean(number_CZ_PCR    (:, 1)) std(number_CZ_PCR    (:, 1)) ...
        mean(number_CNZ_PCR   (:, 1)) std(number_CNZ_PCR   (:, 1));
        mean(number_CZ_ORACLE (:, 1)) std(number_CZ_ORACLE (:, 1)) ...
        mean(number_CNZ_ORACLE(:, 1)) std(number_CNZ_ORACLE(:, 1));
        ];

    RES = [result_RME result_EE result_MSE result_TE result_VS];

    disp(' ');
    disp('  RES = ');
    for i = 1 : size(RES, 1)
        fprintf('%7.2f  %7.4f (%3.1f)  %7.5f   %6.3f   %6.1f (%3.1f)  %4.1f (%3.1f) \n', RES(i, :))
    end
    disp(' ');

    figure(1)
    subplot(2, 2, cases)
    label_vector = (0:1:j_to_display-1)';
    boxplot(bias_beta_PCR(:, :, 1), 'labels', label_vector);
    hold on
    xl = xlim;
    plot([xl(1) xl(2)], [0 0], 'k--'); hold on;
    xlabel('\boldmath$j$', 'Interpreter', 'Latex'); ylabel('');
    if     study == 1
        word_0 = '\textbf{raw}'; %  data

    elseif study == 2
        word_0 = '\textbf{contam.}'; %  data
    end
    if     index_robust(1) == 0 && index_robust(2) == 0
        word_1 = '\textbf{non-robust}';

    elseif index_robust(1) == 1 && index_robust(2) == 1
        word_1 = '\textbf{robust}';
    end
    word_2 = ['{\boldmath$p_n = ', num2str(p_n_set(cases)),'$}'];
    title([word_0, ', ', word_1, ', ', word_2], 'Interpreter', 'Latex');
    ylim([-3 2])

    fprintf('Elapsed time is %f seconds. \n\n', t)
end % for cases = 1:length(num_obs_set)

return;
