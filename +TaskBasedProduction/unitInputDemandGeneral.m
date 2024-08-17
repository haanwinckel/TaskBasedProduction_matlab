function labor_input = unitInputDemandGeneral(xT, q, z, b_g, e_h)
    % unitInputDemandGeneral calculates unit labor demands given an array xT of H-1 thresholds in task space,
    % a productivity value z, a density function b_g for the task distribution, and an array e_h of H functions
    % representing the cost of each labor type as a function of task complexity.
    %
    % Inputs:
    %   xT - A vector of H-1 thresholds in task space (vector)
    %   q - Total production (scalar)
    %   z - Productivity value (scalar)
    %   b_g - A density function for the task distribution (function handle)
    %   e_h - A cell array of H functions representing the cost of each labor type (cell array of function handles)
    %
    % Output:
    %   labor_input - A vector representing the labor demand for each labor type (vector)

    % Ensure xT is a column vector
    xT = xT(:);
    
    % Number of labor types
    H = length(xT) + 1;
    
    % Adding 0 and infinity to thresholds
    xT = [0; xT; Inf];
    
    % Initialize the output vector
    labor_input = zeros(H, 1);
    
    % Check if b_g is a valid density function
    if ~TaskBasedProduction.is_density_function(b_g, xT(1), xT(end))
        error('b_g(x) is not a valid density function');
    end
    
    % Compute labor demand for each type
    for h = 1:H
        integrand = @(x) b_g(x) / (z * e_h{h}(x));
        labor_input(h) = integral(integrand, xT(h), xT(h+1), 'ArrayValued', true);
    end
    
    % Multiply by total production q
    labor_input = q * labor_input;
end
