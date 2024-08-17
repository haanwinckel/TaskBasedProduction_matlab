function [epsilon_h_sub, epsilon_h_compl] = elasticitySubCompGeneral(labor_input, z, b_g, e_h, MPL, xT, q)
    % elasticitySubCompGeneral calculates the elasticity of substitution and complementarity for a given set of parameters.
    %
    % Inputs:
    %   labor_input - Array of labor inputs of different types with H elements. If empty, it will be computed internally given q and xT.
    %   z - Productivity parameter
    %   b_g - General task density function (function handle)
    %   e_h - Cell array of comparative advantage functions (cell array of function handles)
    %   MPL - (optional) Array of marginal products of labor for each worker type with H elements
    %   xT - (optional) Array of task thresholds with H-1 elements
    %   q - (optional) Scalar representing total production. If not provided, it will be computed within the function.
    %
    % Outputs:
    %   epsilon_h_sub - Matrix of elasticity of substitution values (H x H)
    %   epsilon_h_compl - Matrix of elasticity of complementarity values (H x H)

    if nargin < 5 || isempty(MPL), MPL = []; end
    if nargin < 6 || isempty(xT), xT = []; end
    if nargin < 7 || isempty(q), q = []; end

    % If xT or q is empty, compute them using prodFunGeneral
    if isempty(xT) || isempty(q)
        [q, xT] = TaskBasedProduction.prodFunGeneral(labor_input, z, b_g, e_h);
    end

    % If MPL is not provided, compute it
    if isempty(MPL)
        MPL =  TaskBasedProduction.margProdLaborGeneral(labor_input, z, b_g, e_h, xT, q);
    end

    % If labor_input is not provided, compute it using unitInputDemandGeneral
    if isempty(labor_input)
        labor_input =  TaskBasedProduction.unitInputDemandGeneral(xT, q, z, b_g, e_h);
    end

    H = length(labor_input);
    rho_h = zeros(H, 1);
    s_h = (MPL .* labor_input) / q;
    epsilon_h_sub = zeros(H, H);
    epsilon_h_compl = zeros(H, H);
    xT = [xT; Inf];  % Add highest thresholds for the highest worker type

    % Compute effective task outputs and normalized b_g
    e_h_T = arrayfun(@(i) e_h{i}(xT(i)), 1:H);
    b_g_T = arrayfun(@(i) b_g(xT(i)) / z, 1:H);

    % Compute rho_h
    for h = 1:H-1
        % Define the log expression as a function handle
        log_expr = @(x) log(e_h{h+1}(x) / e_h{h}(x));
        
        % Compute the numerical derivative of the log expression at xT(h)
        log_derivative =  TaskBasedProduction.numerical_derivative(log_expr, xT(h));
        
        % Compute rho_h(h)
        rho_h(h) = b_g_T(h) * MPL(h) * (1 / e_h_T(h)) * (1 / log_derivative);
    end

    % Compute epsilon_h_sub
    for h = 1:H
        for h_prime = h+1:H
            if h_prime == h + 1
                epsilon_h_sub(h, h_prime) = rho_h(h) / (s_h(h) * s_h(h_prime));
            end
        end
    end

    % Compute epsilon_h_compl
    xi = zeros(H, H, H);
    temp = zeros(H, H, H);

    for h = 1:H
        for h_prime = h+1:H
            for h_bold = 1:H
                xi(h, h_prime, h_bold) = ((h >= h_bold + 1) - sum(s_h(h_bold + 1:H))) * ((h_bold >= h_prime) - sum(s_h(1:h_bold)));
                temp(h, h_prime, h_bold) = xi(h, h_prime, h_bold) * (1 / rho_h(h_bold));
            end
            epsilon_h_compl(h, h_prime) = sum(temp(h, h_prime, 1:H-1));
        end
    end
end


