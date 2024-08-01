function initial_guess = find_initial_guess(theta, kappa, z, alphaVec, threshold)
    % find_initial_guess generates an initial guess for the optimization problem in prod_fun such that the implied labor demand is non-trivial.
    % Arguments:
    %   theta - The scale parameter of the gamma distribution
    %   kappa - The shape parameter of the gamma distribution
    %   z - A scaling factor for the labor input
    %   alphaVec - An array of task-specific parameters
    %   threshold - The minimum acceptable labor demand for each task (default = 1e-2)
    % Returns:
    %   initial_guess - A vector containing the initial guess for the optimization, including the log of the initial production quantity q and the initial task thresholds xT

    if nargin < 5
        threshold = 1e-2;
    end

    H = length(alphaVec);  % Number of tasks
    initial_q = 0.0;  % log(1) is 0

    % Generate initial xT values using random percentiles from the gamma distribution
    xT = generate_initial_xT(H, kappa, theta);

    % Adjust the xT values to ensure the implied labor demand is non-trivial
    xT = adjust_xT(xT, initial_q, H, theta, kappa, z, alphaVec, threshold);

    % Combine q and xT into the initial guess vector
    initial_guess = [initial_q; xT];
end

function xT = generate_initial_xT(H, kappa, theta)
    % Generate initial xT values using random percentiles from the gamma distribution
    pd = makedist('Gamma', 'a', kappa, 'b', theta);
    percentiles = sort(rand(H-1, 1));  % Generate random percentiles and sort them
    xT = arrayfun(@(p) icdf(pd, p), percentiles);  % Calculate xT for each percentile
end

function xT = adjust_xT(xT, initial_q, H, theta, kappa, z, alphaVec, threshold)
    imp_l = zeros(H, 1);  % Placeholder for labor demand
    max_iterations = 1000;
    iteration = 0;

    while any(imp_l < threshold) && iteration < max_iterations
        try
            imp_xT = cumsum(exp(xT));
            imp_q = exp(initial_q);
            imp_l = TaskBasedProduction.unitInputDemand(imp_xT, imp_q, theta, kappa, z, alphaVec);
        catch
            % If there's an error, generate new initial xT values from scratch
            xT = generate_initial_xT(H, kappa, theta);
            continue;
        end

        if any(imp_l < threshold)
            xT = generate_initial_xT(H, kappa, theta);
        end

        iteration = iteration + 1;
    end

    if iteration == max_iterations
        error('find_initial_guess: Could not find suitable initial xT within the maximum iterations. Consider changing tolerance.');
    end
end
