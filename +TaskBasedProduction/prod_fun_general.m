function [q, xT, fval, initial_guess] = prod_fun_general(labor_input, z, b_g, e_h, varargin)
    % prod_fun_general Calculates the quantity produced and task thresholds
    %
    % Calculates the quantity produced (q), and task thresholds (xT)
    % given labor inputs (labor_input), productivity z, general blueprint density function (b_g),
    % and a vector of efficiency functions (e_h), one for each labor type.
    %
    % Inputs:
    %   labor_input - array of labor inputs of different types (vector)
    %   z - productivity parameter (scalar)
    %   b_g - blueprint density function (function handle)
    %   e_h - vector of efficiency functions, one for each type (cell array of function handles)
    %
    % Optional Inputs:
    %   initial_guess - Initial guess for optimization. If not provided, defaults to zeros array.
    %   x_tol - Tolerance for the solution vector. Default is 1e-12.
    %   f_tol - Tolerance for the function value. Default is 1e-12.
    %   g_tol - Tolerance for the gradient. Default is 1e-12.
    %   iterations - Maximum number of iterations for the optimization. Default is 1000.
    %   max_retries - Maximum number of retries if the optimization fails. Default is 5.
    %   display - Display option for the optimizer. Default is 'off'.
    %   verbose - If true, display detailed information during optimization. Default is false.
    %
    % Returns:
    %   q - quantity produced (scalar)
    %   xT - array of task thresholds (vector)
    %   fval - Final value of the objective function
    %   initial_guess- Initial guess that worked

    % Parse optional inputs
    p = inputParser;
   addParameter(p, 'initial_guess', []);
    addParameter(p, 'x_tol', []);
    addParameter(p, 'f_tol', 1e-6);
    addParameter(p, 'g_tol', []);
    addParameter(p, 'iterations', 1000000);
    addParameter(p, 'max_retries', 1000);
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

    if isempty(initial_guess)
        initial_guess = zeros(length(labor_input) + 1, 1);
    end

    % Set default optimization options for fmincon
    options = optimoptions('fmincon', 'Display', display_option, 'Algorithm', 'interior-point');

    % Override default options if specified
    if ~isempty(x_tol)
        options = optimoptions(options, 'StepTolerance', x_tol);
    end
    if ~isempty(f_tol)
        options = optimoptions(options, 'FunctionTolerance', f_tol);
    end
    if ~isempty(g_tol)
        options = optimoptions(options, 'OptimalityTolerance', g_tol);
    end
    if ~isempty(iterations)
        options = optimoptions(options, 'MaxIterations', iterations);
    end

    % Objective function for optimization
    function val = objFun(x)
        imp_q = exp(x(1));
        imp_xT = cumsum(exp(x(2:end)));
        imp_l = TaskBasedProduction.unitInputDemand_general(imp_xT, imp_q, z, b_g, e_h);
        err = log(imp_l ./ labor_input);
        val = sum(abs(err));
    end

    retry_count = 0;
    success = false;

    % Constraint bounds (modify if you have specific bounds)
    lb = -inf(size(initial_guess)); % Lower bounds
    ub = inf(size(initial_guess));  % Upper bounds

    % Linear constraints (modify if necessary)
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    while retry_count < max_retries && ~success
        try
            % Perform optimization using fmincon
            [x_opt, fval, exitflag] = fmincon(@objFun, initial_guess, A, b, Aeq, beq, lb, ub, [], options);
            if isempty(f_tol)
                f_tol = 1e-4; % Default tolerance
            end
            
            if exitflag > 0 && fval <= f_tol
                q = exp(x_opt(1));
                xT = cumsum(exp(x_opt(2:end)));
                success = true;
                if verbose
                    disp('Optimal Solution found');
                end
            else
                error('prod_fun_general: could not find optimal allocation.');
            end
        catch
            retry_count = retry_count + 1;
            if verbose
                disp('Exit flag:');
                disp(exitflag);
                disp('Function value at solution:');
                disp(fval);
                disp('Function tol:');
                disp(f_tol);
                fprintf('prod_fun_general: An error occurred or could not find optimal allocation. Retrying with a new initial guess. Retry count: %d\n', retry_count);
            end
            initial_guess = TaskBasedProduction.find_initial_guess_gen(z, b_g, e_h, 'threshold', 1e-2);  % Adjust this to your actual method of finding a new initial guess
        end
    end

    if ~success
        error('prod_fun_general: could not find optimal allocation after %d retries.', max_retries);
    end
end