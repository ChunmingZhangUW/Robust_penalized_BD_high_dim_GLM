
function f_lam = true_penalty(x, I_penalty, lambda, deri, options)

%--------------------------------------------------------------------------
% Name    : true_penalty.m
% Function: compute values and derivatives of the true penalty,
%           P_\lambda(|x|), at points they are well-defined
% Penalty : L_1, SCAD, hard thresholding, truncated L_1, MCP
%--------------------------------------------------------------------------
% <Input>
%     x     : argument
% I_penalty : 1 (for L_1), 2 (for SCAD),
%             3 (for hard thresholding), 4 (truncated L_1),
%             5 (for MCP),
%  lambda   : tuning parameter
%   deri    : 0, 1st, 2nd derivative of the penalty
%--------------------------------------------------------------------------
% <Output>
%   f_lam   : values and derivatives of F evaluated at x
%--------------------------------------------------------------------------

f_lam = zeros(size(x)); % initialize

abs_x = abs(x);        % |x|
sig_x = sign(x);       % sign(x)

a = options.SCAD_a;    % 3.7;
P_0 = 0; P_1 = lambda; P_a = a*lambda;

if     I_penalty == 1 % L_1 penalty
    if     deri == 0
        f = lambda*abs_x;

    elseif deri == 1  % discontinuous at P_0
        if abs_x == P_0
            disp('non-exist, return!!!'); return;
        else
            f = lambda*sig_x.*(P_0 < abs_x);
        end

    elseif deri == 2  % discontinuous at P_0
        if abs_x == P_0
            disp('non-exist, return!!!'); return;
        else
            f = 0.*(P_0 < abs_x);
        end
    end
    f_lam = f;

elseif I_penalty == 2 % SCAD penalty
    if     deri == 0
        f_1 = lambda*abs_x.*(abs_x <= P_1);
        f_2 = -((abs_x-P_a).^2-(a^2-1)*lambda^2)...
            /(2*(a-1)).*(P_1 < abs_x & abs_x <= P_a);
        f_3 = (a+1)*lambda^2/2.*(abs_x > P_a);

    elseif deri == 1  % discontinuous at P_0
        if abs_x == P_0
            disp('non-exist, return!!!'); return;
        else
            f_1 = lambda*sig_x.*(P_0 < abs_x & abs_x <= P_1);
            f_2 = -(abs_x-P_a).*sig_x/(a-1).*(P_1 < abs_x & abs_x <= P_a);
            f_3 = 0.*(abs_x > P_a);
        end

    elseif deri == 2  % discontinuous at P_0, P_1 and P_a
        if abs_x == P_0 || abs_x == P_1 || abs_x == P_a
            disp('non-exist, return!!!'); return;
        else
            f_1 = 0.*(P_0 < abs_x & abs_x < P_1);
            f_2 = -1/(a-1).*(P_1 < abs_x & abs_x < P_a);
            f_3 = 0.*(abs_x > P_a);
        end
    end
    f_lam = f_1 + f_2 + f_3;

elseif I_penalty == 3 % hard thresholding penalty
    if     deri == 0
        f_1 = ( lambda^2-(abs_x-P_1).^2 ).*(abs_x <= P_1);
        f_2 = lambda^2.*(abs_x > lambda);

    elseif deri == 1  % discontinuous at P_0
        if abs_x == P_0
            disp('non-exist, return!!!'); return;
        else
            f_1 = -2*(abs_x-P_1).*sig_x.*(P_0 < abs_x && abs_x <= P_1);
            f_2 = 0.*(abs_x > lambda);
        end

    elseif deri == 2  % discontinuous at P_0 and P_1
        if abs_x == P_0 || abs_x == P_1
            disp('non-exist, return!!!'); return;
        else
            f_1 = -2*(P_0 < abs_x & abs_x <= P_1);
            f_2 = 0.*(abs_x > lambda);
        end
    end
    f_lam = f_1 + f_2;

elseif I_penalty == 4 % truncated L_1 penalty
    if     deri == 0
        f_1 = lambda*abs_x.*(abs_x <= P_1);
        f_2 = lambda^2.*(abs_x > P_1);

    elseif deri == 1  % discontinuous at P_0 and P_1
        if abs_x == P_0 || abs_x == P_1
            disp('non-exist, return!!!'); return;
        else
            f_1 = lambda*sig_x.*(P_0 < abs_x & abs_x < P_1);
            f_2 = 0.*(abs_x > P_1);
        end

    elseif deri == 2  % discontinuous at P_0 and P_1
        if abs_x == P_0 || abs_x == P_1
            disp('non-exist, return!!!'); return;
        else
            f_1 = 0.*(P_0 < abs_x & abs_x < P_1);
            f_2 = 0.*(abs_x > P_1);
        end

    end
    f_lam = f_1 + f_2;

elseif I_penalty == 5 % MCP penalty
    if deri == 1  % discontinuous at P_0
        if abs_x == P_0
            disp('non-exist, return!!!'); return;
        else
            f_1 = -(abs_x-P_a).*sig_x/a.*(P_0 < abs_x & abs_x <= P_a);
            f_2 = 0.*(abs_x > P_a);
        end
        f_lam = f_1 + f_2;
    end
end

if any(isnan(f_lam)) == 1
    disp(' !!!true_penalty.m: some estimate of f_lam = NaN!!!');
end

if any(isinf(f_lam)) == 1
    disp(' !!!true_penalty.m: some estimate of f_lam = Inf!!!');
end