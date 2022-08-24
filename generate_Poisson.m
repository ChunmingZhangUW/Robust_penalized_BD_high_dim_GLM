
function [design_matrix, y_vector] = ...
    generate_Poisson(num_obs, p_n, true_beta, family, link, ...
    auxi_set, rho, study, example, A)

%
% Name    : generate_Poisson.m
% Function: generate raw and contaminated data in the simulation study
% Model   : GLM
%--------------------------------------------------------------------------
% <input>
%    num_obs     : # of observations
%    p_n	   : # of coefficients for non-constant variables
% true_beta: (p_n+1)*1 vector, (beta_0,...,beta_d)'
%  family  : 0 (Gaussian), 1 (Bernoulli), 21 (Poisson_quasi)
% auxi_set : set of auxilary parameters
%  example : 1, 2, 3, ... example index
%     A    : 50 or 100, lower bound in the contaminated responses
% <output>
% design_matrix: num_obs*(p_n+1) design matrix, including the column vector (1,...,)'.
% y_vector : num_obs*1 response vector

%Sigma = zeros(p_n-1, p_n-1);

if family == 21 % Poisson_quasi
    mu_vec = zeros(num_obs, p_n-1); % mean vector
    
    cov = rho.^(0:(p_n-2));
    Sigma = toeplitz(cov); % (p_n-1)*(p_n-1) covariance matrix
    %for j = 1:p_n-1
    %    for k = 1:p_n-1
    %        Sigma(j, k) = rho^(abs(j-k)); % (p_n-1)*(p_n-1) covariance matrix
    %    end
    %end
    
    X_d = normcdf(mvnrnd(mu_vec, Sigma));    % map to [0, 1]
    matrix_covariates = [(1:num_obs)'/num_obs-0.5 X_d-0.5];
    % num_obs*p_n matrix
    design_matrix = [ones(num_obs,1) matrix_covariates];
    
    %     mu_vec = zeros(num_obs, p_n); % mean vector
    %     for j = 1:p_n
    %         for k = 1:p_n
    %             Sigma(j, k) = rho^(abs(j-k)); % (p_n-1)*(p_n-1) covariance matrix
    %         end
    %     end
    %     X_d = mvnrnd(mu_vec, Sigma);
    %     X_d = normcdf(X_d);    % map to [0, 1]
    %     design_matrix = [ones(num_obs,1) X_d-0.5];
    
    true_theta = design_matrix * true_beta;
    if strcmpi(link, 'log')
        true_m = exp(true_theta);
    end
    
    true_phi = auxi_set(1);
    if     true_phi > 1   % overdispersed Poisson
        p_vec = 1/true_phi*ones(num_obs,1);         % p in NB(r,p)
        r_vec = true_m/(true_phi-1);          % r in NB(r,p)
        
        ss = 0;
        while ss == 0   % to avoid degenerate cases
            y_vector = nbinrnd(r_vec, p_vec, num_obs, 1);
            ss = sum(y_vector);
        end
        
    elseif true_phi == 1  % classical Poisson
        ss = 0;
        while ss == 0   % to avoid degenerate cases
            y_vector = poissrnd(true_m, num_obs, 1);
            ss = sum(y_vector);
        end
        
    end
    
    if     example == 1    % new
        %------------ contamination -----------------------------
        if     study == 1
            return
            
        elseif study == 2
            U = rand(num_obs, 1); a = 1/2;
            
            design_matrix(1,  2) =  a * sign(U(1)-0.5);
            %y_vector(1) =  10*sign(U(1)-0.5);
            design_matrix(2,  3) =  a * sign(U(2)-0.5);
            %y_vector(2) = -10*sign(U(2)-0.5);
            design_matrix(3,  4) =  a * sign(U(3)-0.5);
            %y_vector(3) =  10*sign(U(3)-0.5);
            design_matrix(4,  6) =  a * sign(U(4)-0.5);
            %y_vector(4) = -10*sign(U(4)-0.5);
            design_matrix(5,  8) =  a * sign(U(5)-0.5);
            %y_vector(5) = -10*sign(U(5)-0.5);
            design_matrix(6,  9) =  a * sign(U(6)-0.5);
            %y_vector(6) = -10*sign(U(6)-0.5);
            design_matrix(7, 10) =  a * sign(U(7)-0.5);
            %y_vector(7) = -10*sign(U(7)-0.5);
            
            %y_vector(1:10) = y_vector(1:10).*(y_vector(1:10) > 100) + 500/5.*(y_vector(1:10) <= 100);
            %y_vector(1:8) = y_vector(1:8).*(y_vector(1:8) > 100) + 100/2.*(y_vector(1:8) <= 100);
            % %version of 11/12/2012
            
            y_vector(1:8) = ...
                y_vector(1:8) .* (y_vector(1:8) >  100) + ...
                A             .* (y_vector(1:8) <= 100);
            % version in 07/??/2013
        else
            disp(' wrong, stop'); return
        end
        %------------------------------------------------------------
        
    elseif example == 2    % old
        %------------ contamination -----------------------------
        if     study == 1
            return
            
        elseif study == 2
            U = rand(num_obs, 1);
            
            a = 3;
            
            design_matrix(1,  2) =  a * sign(U(1)-0.5);
            %y_vector(1) =  10*sign(U(1)-0.5);
            design_matrix(2,  3) =  a * sign(U(2)-0.5);
            %y_vector(2) = -10*sign(U(2)-0.5);
            design_matrix(3,  4) =  a * sign(U(3)-0.5);
            %y_vector(3) =  10*sign(U(3)-0.5);
            design_matrix(4,  6) =  a * sign(U(4)-0.5);
            %y_vector(4) = -10*sign(U(4)-0.5);
            design_matrix(5,  8) =  a * sign(U(5)-0.5);
            %y_vector(5) = -10*sign(U(5)-0.5);
            design_matrix(6,  9) =  a * sign(U(6)-0.5);
            %y_vector(6) = -10*sign(U(6)-0.5);
            design_matrix(7, 10) =  a * sign(U(7)-0.5);
            %y_vector(7) = -10*sign(U(7)-0.5);
            
            y_vector(1:10) = ...
                y_vector(1:10) .* (y_vector(1:10) >  100) + ...
                500            .* (y_vector(1:10) <= 100);
            %y_vector(1:8) = y_vector(1:8).*(y_vector(1:8) > 100) + 100/2.*(y_vector(1:8) <= 100);
            
        else
            disp(' wrong, stop'); return
        end
        %------------------------------------------------------------
        
    end
end

