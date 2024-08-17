function mpl = margProdLaborGeneral(labor_input, z, b_g, e_h, xT, q)
    % margProdLaborGeneral Calculates the marginal productivity of labor for each worker type
    %
    % Calculates the marginal productivity of labor for each worker type given the input parameters.
    %
    % Inputs:
    %   labor_input - Array of labor inputs of different types (vector). If empty, it will be computed internally.
    %   z - Productivity scalar (scalar)
    %   b_g - Task density function (function handle)
    %   e_h - Vector of comparative advantage functions (cell array of function handles)
    %   xT - (optional) Array representing the precomputed task thresholds (vector)
    %   q - (optional) Scalar representing the precomputed quantity produced (scalar)
    %   initial_guess - (optional) Initial guess for the task thresholds, used when computing xT and q.
    %
    % Output:
    %   mpl - Array representing the marginal productivity of labor for each worker type (vector)

    % If xT or q are missing, compute them using prodFunGeneral
    if nargin < 5 || isempty(xT) || isempty(q)
        
        [q, xT] = TaskBasedProduction.prodFunGeneral(labor_input, z, b_g, e_h);
    end

    % If labor_input is missing, calculate it using unitInputDemandGeneral
    if isempty(labor_input)
        labor_input = q * TaskBasedProduction.unitInputDemandGeneral(xT, q, z, b_g, e_h);
    end

    % Number of labor types
    H = length(e_h);
    
    % Pre-allocate array for ratio values
    temp = zeros(H-1, 1); 

    % Calculate the ratio e_{h} / e_{h-1} for h = 2:H and evaluate at xT(h-1)
    for h = 2:H
        ratio_value = e_h{h}(xT(h-1)) / e_h{h-1}(xT(h-1));
        temp(h-1) = ratio_value;
    end

    % Calculate cumulative product for mpl_over_mpl1
    mpl_over_mpl1 = [1; cumprod(temp)];

    % Calculate mpl1
    mpl1 = q / sum(mpl_over_mpl1 .* labor_input);

    % Calculate mpl
    mpl = mpl_over_mpl1 * mpl1;
end


