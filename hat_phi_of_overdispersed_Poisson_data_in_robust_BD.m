
function hat_Poisson_phi = hat_phi_of_overdispersed_Poisson_data_in_robust_BD(...
    I_loss, family, link, theta, y, options, c_tune_constant)

%--------------------------------------------------------------------------
% Name    : hat_phi_of_overdispersed_Poisson_data_in_robust_BD.m
% Function: estimate phi of overdispersed Poisson responses
%           under robust-BD
%--------------------------------------------------------------------------

if I_loss == 211 && family == 21 && strcmpi(link, 'log')
    % V(x)=phi x, (negative) quasi-likelihood; Poisson responses, log link

    mu = exp(theta);

    residual = y - mu;

    %---------------------------------------------------------------------
    Poisson_a_UB = min(( (10^15-mu)/c_tune_constant ).^2./mu);
    square_residual = residual.^2;

    if     options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n == 1

        hat_Poisson_phi = min([...
            mean(square_residual.*mu)/mean(mu.^2), ...
            Poisson_a_UB, options.Poisson_variance_a_UB]);

    elseif options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n == 2

        hat_Poisson_phi = min([max([...
            median(square_residual./mu), ...
            mean(square_residual./mu) ...
            ]), ...
            Poisson_a_UB, options.Poisson_variance_a_UB]);

    elseif options.choice_of_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n == 3
        % specified
        % applicable to the genuine Poisson responses with phi = 1

        %disp(' --> enter hat_phi_of_overdispersed_Poisson_data_in_robust_BD.m: for Poisson with specified phi')
        hat_Poisson_phi = options.specified_Poisson_phi_in_robust_p_1_p_2_BD_for_hat_V_n;

    else
        disp('!!! hat_phi_of_overdispersed_Poisson_data_in_robust_BD.m: hat_Poisson_phi is undefined.')
        return
    end
    %---------------------------------------------------------------------
end
