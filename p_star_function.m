
function f = p_star_function(I_loss, family, theta, y, mu, results, ...
    options, index_robust_c, choice_rho_function, c_tune_constant)

%--------------------------------------------------------------------------
% Name    : p_star_function.m
% Function: explicit formula for p^*(y, mu)
%           = \int_y^mu \psi(r(y,s)) \{q''(s) \sqrt(V(s))\} ds
% Called  : p_star_function_indefinite.m
%--------------------------------------------------------------------------
% <Input>
% index_robust_c:  0: \psi(r)  =  r;
%                  1: \psi(r) \ne r;
%choice_rho_function: 1 Huber \rho function; 2 Tukey biweight function
% c_tune_constant: constant used in \psi(r)
%--------------------------------------------------------------------------
% <Output>
%--------------------------------------------------------------------------

%--------------------- for classical-BD ----------------------

if index_robust_c == 0 || c_tune_constant == inf  % \psi(r) = r

    f = true_qq(I_loss, family, link, theta, y, deri, options);
    return
end

%---------------------- for robust-BD ------------------------

% index_robust_c == 1 && c_tune_constant ~= inf  % \psi(r) \ne r

p_star_y_mu = p_star_function_indefinite(I_loss, family, theta, y, mu, ...
    results, options, index_robust_c, choice_rho_function, c_tune_constant);

if     family == 12  % Binomial responses
    y_Bin = y(:, 1);
    p_star_y_y  = p_star_function_indefinite(I_loss, family, theta, y, y_Bin, ...
        results, options, index_robust_c, choice_rho_function, c_tune_constant);

elseif family ~= 12  % non-Binomial responses
    p_star_y_y  = p_star_function_indefinite(I_loss, family, theta, y, y, ...
        results, options, index_robust_c, choice_rho_function, c_tune_constant);
end

f = p_star_y_mu - p_star_y_y;

if any(isempty(f)) == 1
    disp(' !!!p_star_function.m: some estimate of f = [ ]!!!');
end

if any(isnan(f)) == 1
    disp(' !!!p_star_function.m: some estimate of f = NaN!!!');
end

if any(isinf(f)) == 1
    disp(' !!!p_star_function.m: some estimate of f = Inf!!!');
end
