
function [q_2_vector, q_1_d_q_2] = weight_response_BD(...
    I_loss, family, link, theta, y, options)

%--------------------------------------------------------------------------
% Name    : weight_response_BD.m
% Function: compute qq_2(y;\theta) and qq_1(y;\theta)/qq_2(y;\theta)
%           used in weighted quadratic approximation
%           to Bregman Divergence Q_q(y, F^{-1}(\theta))
% Model   : F(m(x))=beta_1 x_1+...+beta_K x_K, where m(x)=E(Y|X=x)
% Loss    : deviance loss, exponential loss, quadratic loss,
%           quasi-likelihood
% Link F  : the canonical link function
% Used in : weighted quadratic approximation
%--------------------------------------------------------------------------
% <Input>
% I_loss : choice of loss function:
%          1 (deviance), 2 (exponential), 3 (quadratic), 4 (arching),
%          211 (V(x)=phi x)
% family : 0 (Gaussian); 1 (Bernoulli), 12 (Binomial); 21 (Poisson_quasi)
%  link  : type of link function, 'iden', 'logit', 'log'
%  theta : scalar or vector
%    y   : scalar or vector with the same size of theta
% options: set of parameters
%--------------------------------------------------------------------------
% <Output>
% q_2_vector: same size of y
% q_1_d_q_2 : same size of y
%--------------------------------------------------------------------------

if     family == 0 && strcmpi(link, 'iden')
    % Gaussian responses, identity link

    mu = theta;

    if (I_loss == 1 || I_loss == 3) % deviance, quadratic
        q_2_vector = 2*ones(length(y), 1);

        q_1_d_q_2 = -(y - mu);
    end

elseif family == 1 && strcmpi(link, 'logit')
    % Bernoulli responses, logit link

    mu = 1./(1 + exp(-theta));

    mu_modified           = max(    mu, options.zero_thres);
    one_minus_mu_modified = max(1 - mu, options.zero_thres);

    if     I_loss == 1 % deviance
        q_2_vector = 2*mu.*(1 - mu);

        ratio = ...
            + 1./mu_modified.*y ...
            - 1./one_minus_mu_modified.*(1 - y);
        q_1_d_q_2 = -ratio;

    elseif I_loss == 2 % exponential
        q_2_vector = ( ...
            + sqrt(mu./one_minus_mu_modified).*(1 - y) ...
            + sqrt((1 - mu)./mu_modified).*y )/4;

        q_1_d_q_2 = 2 - 4*y; % -4*(y-1/2);
    end

elseif family == 12 && strcmpi(link, 'logit')
    % Binomial responses, logit link

    y_Bin = y(:, 1);
    N_Bin = y(:, 2);

    p_Bin = 1./(1 + exp(-theta));

    mu = N_Bin.*p_Bin;

    mu_modified         = max(        mu, options.zero_thres);
    N_minus_mu_modified = max(N_Bin - mu, options.zero_thres);

    V_mu_original = mu.*(1 - p_Bin);    % mu.*(1-p);

    if I_loss == 1 % deviance
        q_2_vector = 2*V_mu_original;

        ratio = ...
            + N_Bin./mu_modified         .*(y_Bin == N_Bin) ...
            - N_Bin./N_minus_mu_modified .*(y_Bin == 0) ...
            + (y_Bin - mu).*N_Bin./max(mu.*(N_Bin - mu), options.zero_thres) ...
            .*(0 < y_Bin & y_Bin < N_Bin);
        q_1_d_q_2 = -ratio;
    end

elseif family == 21 && strcmpi(link, 'log')
    % Poisson_quasi, log link

    theta(theta < -options.BD_C) = -options.BD_C;
    theta(theta >  options.BD_C) =  options.BD_C;

    mu = exp(theta);

    mu_modified = max(mu, options.zero_thres);

    if I_loss == 211  % V(x)=phi x, (negative) quasi-likelihood
        q_2_vector = mu;

        q_1_d_q_2 = 1 - y./mu_modified; % -(y - mu)./mu;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q_2_vector = max(q_2_vector, options.delta);  % set to be >= options.delta

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if any(isempty(q_1_d_q_2)) == 1
    disp(' !!!weight_response_BD.m: some estimate of q_1_d_q_2 = [ ]!!!');
end

if any(isnan(q_1_d_q_2))   == 1
    disp(' !!!weight_response_BD.m: some estimate of q_1_d_q_2 = NaN!!!');
end

if any(isinf(q_1_d_q_2))   == 1
    disp(' !!!weight_response_BD.m: some estimate of q_1_d_q_2 = Inf!!!');
end
%-------------

if any(isempty(q_2_vector)) == 1
    disp(' !!!weight_response_BD.m: some estimate of q_2_vector = [ ]!!!');
end

if any(isnan(q_2_vector))   == 1
    disp(' !!!weight_response_BD.m: some estimate of q_2_vector = NaN!!!');
end

if any(isinf(q_2_vector))   == 1
    disp(' !!!weight_response_BD.m: some estimate of q_2_vector = Inf!!!');
end
