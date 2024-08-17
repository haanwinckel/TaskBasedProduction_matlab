function labor_input = unitInputDemand(xT, q, theta, kappa, z, alphaVec, skipParamChecks) 
    % Calculates unit labor demands given blueprint scale theta, blueprint shape kappa,
    % productivity z, an array of comparative advantage values alphaVec with H elements
    % (one for each worker type), and an array xT of H-1 thresholds in task space.
    %
    % Inputs:
    %   xT - array of thresholds in task space (vector)
    %   q - Total production (scalar)
    %   theta - blueprint scale (scalar)
    %   kappa - blueprint shape (scalar)
    %   z - productivity (scalar)
    %   alphaVec - array of comparative advantage values (vector)
    %   skipParamChecks - flag to skip parameter checks (boolean)
    %
    % Output:
    %   labor_input - unit labor demands (vector)

    % Assign default value to skipParamChecks if not provided
    if nargin < 7
        skipParamChecks = false;
    end

    % Parameter checks if not skipped
    if ~skipParamChecks
        assert(theta > 0, 'theta must be greater than 0');
        assert(kappa > 0, 'kappa must be greater than 0');
        assert(z > 0, 'z must be greater than 0');
        assert(all(diff(alphaVec) > 0), 'alphaVec values must be increasing');
        assert(all(diff(xT) >= 0), 'xT values must be non-decreasing');
    end

    % Number of labor types
    H = length(xT) + 1;
    
    % Initialize labor input as a zero vector
    labor_input = zeros(H, 1);
    
    % Add 0 and Inf to the threshold array xT
    xT = [0.0; xT(:); Inf]; 

    % Loop through each labor type
    for h = 1:H
        alpha = alphaVec(h);  % Get the comparative advantage value for the h-th type
        upsilon = alpha + 1/theta;  % Calculate upsilon

        % Depending on the value of upsilon, calculate labor demand
        if upsilon > 0.0
            % Positive upsilon case
            labor_input(h) = TaskBasedProduction.component_positive_ups(upsilon, kappa, xT(h), xT(h+1)) / upsilon^kappa;
        elseif upsilon == 0.0
            % Zero upsilon case
            labor_input(h) = (xT(h+1)^kappa - xT(h)^kappa) / (kappa * gamma(kappa));
        else
            % Negative upsilon case
            labor_input(h) = TaskBasedProduction.component_negative_ups(upsilon, kappa, xT(h), xT(h+1));
        end
    end

    % Scale labor_input by the production quantity, productivity, and blueprint parameters
    labor_input = q * labor_input / (z * theta^kappa);
end