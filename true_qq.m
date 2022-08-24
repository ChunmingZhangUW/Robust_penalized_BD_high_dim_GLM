
function f = true_qq(I_loss, family, link, theta, y, deri, options)

%--------------------------------------------------------------------------
% Name    : true_qq.m
% Function: compute qq_j(y; \theta), i.e. values and derivatives of
%           Bregman Divergence Q_q(y, F^{-1}(\theta)) w.r.t. \theta
% Loss    : deviance loss, exponential loss, quadratic loss,
%           quasi-likelihood
%--------------------------------------------------------------------------
% <Input>
% I_loss : choice of loss function:
%          1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%          211 (V(x)=phi x)
% family : 0 (Gaussian); 1 (Bernoulli), 12 (Binomial); 21 (Poisson_quasi)
%  link  : type of link function, 'iden', 'logit', 'log'
%  theta : scalar or vector
%    y   : scalar or vector with the same size of theta
%   deri : 0, 1, 2 derivatives of the loss
% options: set of parameters
%--------------------------------------------------------------------------
% <Output>
%    f   : values and derivatives of the BD Q-loss, same size of y
%--------------------------------------------------------------------------

n_obs = size(y, 1);  % y is n_obs*2 for Binomial response dataset

if     family == 0 && strcmpi(link, 'iden')
    % Gaussian responses, identity link

    mu = theta;

    if (I_loss == 1 || I_loss == 3) % deviance, quadratic
        if     deri == 0 % quadratic
            q_0_vector = (y - mu).^2;

            f = q_0_vector;

        elseif deri == 1 % quadratic
            q_1_vector = -2*(y - mu);

            f = q_1_vector;

        elseif deri == 2 % quadratic
            q_2_vector = 2*ones(n_obs, 1);

            f = q_2_vector;
        end
    end

elseif family == 1 && strcmpi(link, 'logit')
    % Bernoulli responses, logit link

    mu = 1./(1 + exp(-theta));

    if     I_loss == 1 % deviance
        if     deri == 0 % deviance
            mu_modified           = max(    mu, options.zero_thres);
            one_minus_mu_modified = max(1 - mu, options.zero_thres);

            log_mu_modified           = log(mu_modified);
            log_one_minus_mu_modified = log(one_minus_mu_modified);

            q_0_vector = ...
                -2*log_one_minus_mu_modified.*(1 - y) ...
                -2*log_mu_modified          .*y;  % more stable than the one below
            %q_0_vector = -2*(y.*theta + log_one_minus_mu_modified);

            f = q_0_vector;

        elseif deri == 1 % deviance
            q_1_vector = -2*(y - mu);

            f = q_1_vector;

        elseif deri == 2 % deviance
            if     options.UB_Bernoulli_dev_choice == 0
                q_2_vector = 2*mu.*(1 - mu);

            elseif options.UB_Bernoulli_dev_choice == 1
                q_2_vector = ones(n_obs, 1)/2;
            end

            f = q_2_vector;
        end

    elseif I_loss == 2 % exponential
        mu_modified           = max(    mu, options.zero_thres);
        one_minus_mu_modified = max(1 - mu, options.zero_thres);

        if     deri == 0 % exponential
            q_0_vector = ...
                + sqrt(mu./one_minus_mu_modified).*(1 - y) ...
                + sqrt((1 - mu)./mu_modified)    .*y;

            f = q_0_vector;

        elseif deri == 1 % exponential
            ratio = ...
                - sqrt(mu./one_minus_mu_modified).*(1 - y) ...
                + sqrt((1 - mu)./mu_modified)    .*y;
            q_1_vector = -ratio/2;

            f = q_1_vector;

        elseif deri == 2 % exponential
            q_2_vector = ( ...
                + sqrt(mu./one_minus_mu_modified).*(1 - y) ...
                + sqrt((1 - mu)./mu_modified)    .*y ...
                )/4;

            f = q_2_vector;
        end
    end

elseif family == 12 && strcmpi(link, 'logit')
    % Binomial responses, logit link

    y_Bin = y(:, 1);
    N_Bin = y(:, 2);

    p_Bin = 1./(1 + exp(-theta));

    mu = N_Bin.*p_Bin;                % N_Bin*p_Bin

    V_mu_original = mu.*(1 - p_Bin);  % N_Bin*p_Bin*(1-P_Bin)= mu.*(1-p);

    if I_loss == 1 % deviance
        if     deri == 0
            N_minus_y  = N_Bin - y_Bin;    % N - y
            N_minus_mu = N_Bin - mu;       % N - mu

            mu_modified         = max(        mu, options.zero_thres);
            N_minus_mu_modified = max(N_minus_mu, options.zero_thres);

            I_2 = -2*( ...
                + y_Bin    .*log(mu_modified) ...
                + N_minus_y.*log(N_minus_mu_modified) );
            % -2\{y log(u) + (N - y) log(N - u)\}
            %------------------

            N_log_N = N_Bin.*log(N_Bin);   % N log(N)

            I_1 = zeros(n_obs, 1);
            I_1(y_Bin == 0    ) = -2*N_log_N(y_Bin == 0    );
            I_1(y_Bin == N_Bin) = -2*N_log_N(y_Bin == N_Bin);
            I_1(0 < y_Bin & y_Bin < N_Bin) = -2*( ...
                + y_Bin         (0 < y_Bin & y_Bin < N_Bin)   ...
                .*log( y_Bin    (0 < y_Bin & y_Bin < N_Bin) ) ...
                + N_minus_y     (0 < y_Bin & y_Bin < N_Bin)   ...
                .*log( N_minus_y(0 < y_Bin & y_Bin < N_Bin) ) ...
                );
            % -2\{y log(y) + (N - y) log(N - y)\}
            %------------------

            q_0_vector = I_2 - I_1;

            f = q_0_vector;

        elseif deri == 1
            q_1_vector = -2*(y_Bin - mu);

            f = q_1_vector;

        elseif deri == 2
            q_2_vector = 2*V_mu_original;

            f = q_2_vector;
        end
    end

elseif family == 21 && strcmpi(link, 'log')
    % Poisson_quasi, log link

    theta(theta < -options.BD_C) = -options.BD_C;
    theta(theta >  options.BD_C) =  options.BD_C;

    mu = exp(theta);

    if I_loss == 211  % V(x)=phi x, (negative) quasi-likelihood
        if     deri == 0
            y_log_y = zeros(n_obs, 1);  % y*log(y)
            %y_log_y(y == 0) = 0;
            y_log_y(y > 0) = y(y > 0).*log(y(y > 0));

            I_2 = -(y.*theta - mu);
            I_1 = -(y_log_y - y);

            q_0_vector = I_2 - I_1;

            f = q_0_vector;

        elseif deri == 1
            q_1_vector = -(y - mu);

            f = q_1_vector;

        elseif deri == 2
            q_2_vector = mu;

            f = q_2_vector;
        end
    end
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if any(isempty(f)) == 1
    disp([' !!!true_qq.m: some estimate of f = [ ]!!!', ', I_loss=', ...
        num2str(I_loss), ', deri = ', num2str(deri)]);
end

if any(isnan(f)) == 1
    disp([' !!!true_qq.m: some estimate of f = NaN!!!', ', I_loss=', ...
        num2str(I_loss), ', deri = ', num2str(deri)]);
end

if any(isinf(f)) == 1
    disp([' !!!true_qq.m: some estimate of f = Inf!!!', ', I_loss=', ...
        num2str(I_loss), ', deri = ', num2str(deri)]);
end
