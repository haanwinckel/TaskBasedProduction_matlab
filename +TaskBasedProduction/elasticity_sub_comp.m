function [epsilon_h_sub, epsilon_h_compl] = elasticity_sub_comp(labor_input, theta, kappa, z, alphaVec, MPL, xT, initial_guess)
    % elasticity_sub_comp calculates the elasticity of substitution and complementarity for a given set of parameters.
    % Inputs:
    %   labor_input - Array of labor inputs of different types with H elements
    %   theta - Blueprint scale parameter
    %   kappa - Blueprint shape parameter
    %   z - Productivity parameter
    %   alphaVec - Array of comparative advantage values with H elements
    %   MPL - (optional) Array representing the marginal productivity of labor. If not provided, it will be computed within the function.
    %   xT - (optional) Array representing precomputed task thresholds. If not provided, it will be computed within the function.
    %   initial_guess - (optional) Initial guess for the task thresholds, used when computing xT.
    % Outputs:
    %   epsilon_h_sub - Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns)
    %   epsilon_h_compl - Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns)

    if nargin < 6 || isempty(MPL)
        MPL = TaskBasedProduction.margProdLabor(labor_input, theta, kappa, z, alphaVec);  % Compute MPL if not provided
    end
    
    if nargin < 7 || isempty(xT)
        if isempty(initial_guess)
            initial_guess=TaskBasedProduction.find_initial_guess(theta, kappa, z, alphaVec);
        end
        [q, xT, fval] = TaskBasedProduction.prod_fun(labor_input, theta, kappa, z, alphaVec, 'initial_guess', initial_guess);
    end

    H = length(alphaVec);
    rho_h = zeros(H, 1);
    s_h = (MPL .* labor_input) / sum(MPL .* labor_input);
    epsilon_h_sub = zeros(H, H);
    epsilon_h_compl = zeros(H, H);
    xT = [xT; Inf];  % Add highest threshold for highest worker type

    e_h_T = exp(alphaVec .* xT);  % Denominator of rho_h

    for h = 1:H-1
        b_g_T = xT(h)^(kappa - 1) * (1 / (z * gamma(kappa) * theta^kappa)) * exp(-xT(h) / theta);
        rho_h(h) = b_g_T * MPL(h) * (1 / e_h_T(h)) * (1 / (alphaVec(h+1) - alphaVec(h)));
    end

    for h = 1:H
        for h_prime = h+1:H
            if h_prime == h + 1
                epsilon_h_sub(h, h_prime) = rho_h(h) / (s_h(h) * s_h(h_prime));
            end
        end
    end

    xi = zeros(H, H, H);
    temp = zeros(H, H, H);

    for h = 1:H
        for h_prime = h+1:H
            for h_bold = 1:H
                xi(h, h_prime, h_bold) = ((h >= h_bold + 1) - sum(s_h(h_bold + 1:H))) * ((h_bold >= h_prime) - sum(s_h(1:h_bold)));
                temp(h, h_prime, h_bold) = xi(h, h_prime, h_bold) * (1 / rho_h(h_bold));
            end
            if h < h_prime
                epsilon_h_compl(h, h_prime) = sum(temp(h, h_prime, 1:H-1));
            end
        end
    end
end
