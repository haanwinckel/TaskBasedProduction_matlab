function [q, xT] = prodFunGeneral(labor_input, z, b_g, e_h, varargin)
    % prodFunGeneral Calculates the quantity produced and task thresholds
    %
    % Calculates the quantity produced (q) and task thresholds (xT)
    % given labor inputs (labor_input), productivity z, general blueprint density function (b_g),
    % and a vector of efficiency functions (e_h), one for each labor type.
    %
    % Inputs:
    %   labor_input - Array of labor inputs of different types (vector)
    %   z - Productivity parameter (scalar)
    %   b_g - Blueprint density function (function handle)
    %   e_h - Vector of efficiency functions, one for each labor type (cell array of function handles)
    %
    % Optional Inputs:
    %   initial_guess - Initial guess for optimization. If not provided, it will be generated.
    %   x_tol - Tolerance for the solution vector. Default is 1e-6.
    %   f_tol - Tolerance for the function value. Default is 1e-6.
    %   g_tol - Tolerance for the gradient. Default is 1e-6.
    %   iterations - Maximum number of iterations for the optimization. Default is 1000.
    %   max_retries - Maximum number of retries if the optimization fails. Default is 50.
    %   display - Display option for the optimizer. Default is 'off'.
    %   verbose - If true, display detailed information during optimization. Default is false.
    %
    % Outputs:
    %   q - Quantity produced (scalar)
    %   xT - Array of task thresholds (vector)
    %   fval - Final value of the objective function
    %   initial_guess - Initial guess that worked

    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'initial_guess', []);
    addParameter(p, 'x_tol', 1e-6);
    addParameter(p, 'f_tol', 1e-6);
    addParameter(p, 'g_tol', 1e-6);
    addParameter(p, 'iterations', 1000);
    addParameter(p, 'max_retries', 50);
    addParameter(p, 'display', 'off');
    addParameter(p, 'verbose', false);
    parse(p, varargin{:});
    
    initial_guess = p.Results.initial_guess;
    x_tol = p.Results.x_tol;
    f_tol = p.Results.f_tol;
    g_tol = p.Results.g_tol;
    iterations = p.Results.iterations;
    max_retries = p.Results.max_retries;
    display_option = p.Results.display;
    verbose = p.Results.verbose;

    % If no initial guess is provided, generate it using getStartGuessGen_xT
    if isempty(initial_guess)
        initial_guess = TaskBasedProduction.getStartGuessGen_xT(z, b_g, e_h);
    end

    % Set default optimization options for fmincon
    options = optimoptions('fmincon', 'Display', display_option, 'Algorithm', 'interior-point', ...
        'StepTolerance', x_tol, 'FunctionTolerance', f_tol, 'OptimalityTolerance', g_tol, ...
        'MaxIterations', iterations);

    % Objective function for optimization
    function val = objFun(x)
        imp_q = exp(x(1));
        imp_xT = cumsum(exp(x(2:end)));
        imp_l = TaskBasedProduction.unitInputDemandGeneral(imp_xT, imp_q, z, b_g, e_h);
        err = log(imp_l ./ labor_input);
        val = sum(abs(err));  % Objective function is the sum of absolute errors
    end

    retry_count = 0;
    success = false;

    % Constraint bounds (modify if necessary)
    lb = -inf(size(initial_guess));  % Lower bounds
    ub = inf(size(initial_guess));   % Upper bounds

    % Optimization loop with retry logic
    while retry_count < max_retries && ~success
        try
            % Perform optimization using fmincon
            [x_opt, fval, exitflag] = fmincon(@objFun, initial_guess, [], [], [], [], lb, ub, [], options);
            
            % Check if the optimization succeeded
            if exitflag > 0 && fval <= f_tol
                q = exp(x_opt(1));
                xT = cumsum(exp(x_opt(2:end)));
                success = true;
                if verbose
                    disp('prodFunGeneral: Optimal solution found.');
                end
            else
                error('prodFunGeneral: Could not find optimal allocation.');
            end

        catch
            retry_count = retry_count + 1;
            if verbose
                fprintf('prodFunGeneral: Error occurred. Retrying with a new initial guess. Retry count: %d\n', retry_count);
            end
            initial_guess = TaskBasedProduction.getStartGuessGen_xT(z, b_g, e_h);
        end
    end

    if ~success
        error('prodFunGeneral: Could not find optimal allocation after %d retries.', max_retries);
    end
end