
function deri_1_G_func = deri_1_G(family, I_loss, mu, ...
    index_robust_c, choice_rho_function, c_tune_constant, ...
    y, results, options)

%--------------------------------------------------------------------------
% Name    : deri_1_G.m
% Function: compute G'(mu) as a function of mu
% Called  : deri_G_1_function.m
%--------------------------------------------------------------------------
% <Input>
%      mu       : scalar
% index_robust_c:  0: \psi(r)  =  r;
%                  1: \psi(r) \ne r;
%choice_rho_function: 1 for Huber \rho function; 2 for Tukey biweight function
% c_tune_constant: constant used in \psi(r)
%      y        : scalar, or 1\times 2 for the Binomimal responses
%                 y is included to obtain N_Bin for Binomial responses
%--------------------------------------------------------------------------
% <Output>
% deri_1_G_func : scalar
%--------------------------------------------------------------------------

%--------------------- for classical-BD ----------------------

if index_robust_c == 0 || c_tune_constant == inf  % \psi(r) = r

    deri_1_G_func = 0;
    return
end

%---------------------- for robust-BD ------------------------

% index_robust_c == 1 && c_tune_constant ~= inf  % \psi(r) \ne r

N = length(mu);
mu = reshape(mu, N, 1);

if     family == 0     % Gaussian responses
    V_original = results.V_original;
    sqrt_V_original = sqrt(V_original);

    if (I_loss == 1 || I_loss == 3) % deviance, quadratic
        deri_2_q = -2;
    end

elseif family == 1     % Bernoulli responses
    one_minus_mu = 1 - mu;
    V_original = mu.*one_minus_mu;
    sqrt_V_original = sqrt(V_original);
    V_modified      = max(options.zero_thres, V_original);

    if     I_loss == 1  % deviance
        deri_2_q = -2./V_modified;

    elseif I_loss == 2  % exponential
        deri_2_q = -1./(2*V_modified.^(3/2));
    end

elseif family == 12    % Binomial responses
    N_Bin = y(:, 2);

    p_Bin = mu./N_Bin;

    V_original = mu.*(1 - p_Bin);    % mu.*(1-p);
    sqrt_V_original = sqrt(V_original);
    V_modified      = max(options.zero_thres, V_original);

    if I_loss == 1  % deviance
        deri_2_q = -2./V_modified;
    end

elseif family == 21    % Poisson responses
    Poisson_phi = results.Poisson_phi;

    mu_modified = max(options.zero_thres, mu);

    V_original = Poisson_phi*mu;  % V(mu) = \phi mu
    sqrt_V_original = sqrt(V_original);

    if I_loss == 211  % V(x)=phi x, (negative) quasi-likelihood
        deri_2_q = -1./mu_modified;
    end
end

deri_G_1 = deri_G_1_function(family, mu, y, results, options, ...
    index_robust_c, choice_rho_function, c_tune_constant, 1);
% G_1'(mu)

deri_1_G_func = deri_G_1.*(deri_2_q.*sqrt_V_original);

if any(isempty(deri_1_G_func)) == 1
    disp(' !!!deri_1_G.m: some estimate of deri_1_G_func = [ ]!!!');
end

if any(isnan(deri_1_G_func)) == 1
    disp(' !!!deri_1_G.m: some estimate of deri_1_G_func = NaN!!!');
end

if any(isinf(deri_1_G_func)) == 1
    disp(' !!!deri_1_G.m: some estimate of deri_1_G_func = Inf!!!');
end
