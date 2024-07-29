function initial_guess = find_initial_guess_gen(z, b_g, e_h, varargin)
    % find_initial_guess_gen generates an initial guess for the optimization problem
    % using a general density function such that the implied labor demand is non-trivial.
    %
    % Inputs:
    %   z - A scaling factor for the labor input (scalar)
    %   b_g - The general density function (function handle)
    %   e_h - An array of task-specific functions representing the cost of each labor type (cell array of function handles)
    %
    % Optional Inputs:
    %   threshold - The minimum acceptable labor demand for each task (default is 1e-2)
    %   verbose - Boolean flag to enable or disable verbose output for debugging (default is false)
    %
    % Output:
    %   initial_guess - A vector containing the initial guess for the optimization, 
    %                   including the log of the initial production quantity q and the initial task thresholds xT

    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'threshold', 1e-2);
    addParameter(p, 'verbose', false);
    parse(p, varargin{:});
    
    threshold = p.Results.threshold;
    verbose = p.Results.verbose;
    
    H = length(e_h);  % Number of labor types
    initial_q = 0.0;  % log(1) is 0
    
    % Define the CDF by integrating the PDF
    function cdf_val = c_cdf(b_g, x)
        cdf_val = integral(b_g, 0, x);
    end

    % Function to calculate the percentile using the general CDF
    function x_percentile = calculate_percentile(p, b_g)
        % Define the error function for optimization
        function err = err_func(x)
            if verbose
                fprintf('Evaluating CDF at x = %f\n', x);
            end
            if x < 0
                err = Inf;  % Return a large error for invalid domain values
                return;
            end
            cdf_val = c_cdf(b_g, x);
            if verbose
                fprintf('CDF value: %f, Target percentile: %f\n', cdf_val, p);
            end
            err = (cdf_val - p)^2;  % Squared difference to ensure non-negative optimization
        end

        % Perform the optimization using fminunc
        initial_guess = 0.1;  % Initial guess
        options = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton');
        x_percentile = fminunc(@err_func, initial_guess, options);
    end

    % Function to generate initial xT values using random percentiles from the general CDFs
    function xT = generate_initial_xT()
        percentiles = sort(rand(H-1, 1));  % Generate random percentiles and sort them
        xT = arrayfun(@(p) calculate_percentile(p, b_g), percentiles);  % Calculate xT for each percentile
    end

    function xT = adjust_xT(xT)
        imp_l = zeros(H, 1);  % Placeholder for labor demand
        max_iterations = 1000;
        iteration = 0;

        while any(imp_l < threshold) && iteration < max_iterations
            try
                imp_xT = cumsum(exp(xT));
                imp_l = exp(initial_q) * TaskBasedProduction.unitInputDemand_general(imp_xT, z, b_g, e_h);
            catch
                % If there's an error, generate new initial xT values from scratch
                xT = generate_initial_xT();
                continue;
            end

            if any(imp_l < threshold)
                xT = generate_initial_xT();
            end
           
            iteration = iteration + 1;
        end

        if iteration == max_iterations
            error('find_initial_guess_gen: Could not find suitable initial xT within the maximum iterations. Consider changing tolerance.');
        end
    end

    initial_xT = generate_initial_xT();
    adjusted_xT = adjust_xT(initial_xT);

    % Combine q and xT into the initial guess vector
    initial_guess = [initial_q; adjusted_xT(:)];
end