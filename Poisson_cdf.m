
function F = Poisson_cdf(J_vec, mu_vec)
%--------------------------------------------------------------------------
% Name    : Poisson_cdf.m
% Function: compute the cdf function of the Poisson distribution
% Called  : poisscdf.m
%--------------------------------------------------------------------------

F = poisscdf(J_vec, mu_vec);

return

%-------------------------------------------------
nn = length(J_vec);
F = zeros(nn, 1);      % initialize

for k = 1:nn

    j = J_vec(k);      % P(X=j)
    mu_j = mu_vec(k);
    P_0 = exp(-mu_j);

    if     j == 0
        F(k) = P_0;

    elseif j >= 1
        if j == inf
            F(k) = 1; % F(+\infty)=1

        else
            %             P_L = P_0; F(k) = P_0;
            %
            %             for L = 1 : min(j, UB)
            %                 P_L = P_L * mu_j/L;
            %
            %                 F(k) = F(k) + P_L;
            %             end
            F(k) = poisscdf(j, mu_j);

            if any(isempty(F(k))) == 1
                disp([' !!!Poisson_cdf.m, poisscdf(',num2str(j),') = [ ]!!!']);
            end

            if any(isnan(F(k)))   == 1
                disp([' !!!Poisson_cdf.m, poisscdf(',num2str(j),') = NaN!!!']);
            end

            if any(isinf(F(k)))   == 1
                disp([' !!!Poisson_cdf.m, poisscdf(',num2str(j),') = Inf!!!']);
            end

        end
    end

end

F = min(max(0, F), 1);   % CDF is in [0, 1]

return

