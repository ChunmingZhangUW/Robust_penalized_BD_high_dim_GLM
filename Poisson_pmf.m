
function f = Poisson_pmf(J_vec, mu_vec)
%--------------------------------------------------------------------------
% Name    : Poisson_pmf.m
% Function: compute the pmf function of the Poisson distribution
% Called  : poisspdf.m
%--------------------------------------------------------------------------

f = poisspdf(J_vec, mu_vec);

return

%-------------------------------------------------
nn = length(J_vec);
f = zeros(nn, 1);      % initialize

for k = 1:nn

    j = J_vec(k);      % P(X=j) for X ~ Poisson(mu_k)
    mu_k = mu_vec(k);
    P_0 = exp(-mu_k);

    if     j == 0
        f(k) = P_0;

    elseif j >= 1
        if j == inf
            f(k) = 0; % f(+\infty)=0

        else
            %             P_L = P_0;
            %
            %             for L = 1 : min(j, UB)
            %                 P_L = P_L * mu_k/L;
            %             end
            %
            %             f(k) = P_L;
            f(k) = poisspdf(j, mu_k);

            if any(isempty(f(k))) == 1
                disp([' !!!Poisson_pdf.m, poisspdf(',num2str(j),') = [ ]!!!']);
            end

            if any(isnan(f(k)))   == 1
                disp([' !!!Poisson_pdf.m, poisspdf(',num2str(j),') = NaN!!!']);
            end

            if any(isinf(f(k)))   == 1
                disp([' !!!Poisson_pdf.m, poisspdf(',num2str(j),') = Inf!!!', ...
                    ', set it = 0!!!']);
                f(k) = 0;
            end

        end
    end

end

f = min(max(0, f), 1);   % PMF is in [0, 1]

return
% poisspdf   (floor(6.2956*(1.0e+021)+1), 6.2956*(1.0e+021)) = inf, under
%                 MATLAB Version 7.8.0.347 (R2009a)
% Poisson_pmf(floor(6.2956*(1.0e+021)+1), 6.2956*(1.0e+021)) = 0
