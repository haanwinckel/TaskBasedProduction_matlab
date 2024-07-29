function mpl = margProdLabor_general(labor_input, z, b_g, e_h, xT, q, initial_guess)
    % margProdLabor_general Calculates the marginal productivity of labor for each worker type
    %
    % Calculates the marginal productivity of labor for each worker type given the input parameters.
    %
    % Inputs:
    %   labor_input - array of labor inputs of different types (vector)
    %   z - productivity scalar (scalar)
    %   b_g - task density function (function handle)
    %   e_h - vector of comparative advantage functions (cell array of function handles)
    %   xT - (optional) array representing the precomputed task thresholds (vector)
    %   q - (optional) scalar representing the precomputed quantity produced (scalar)
    %   initial_guess - (optional) Initial guess for the task thresholds, used when computing xT and q.
    %
    % Output:
    %   mpl - array representing the marginal productivity of labor for each worker type (vector)

    if nargin < 5 || isempty(xT) || isempty(q)
        if isempty(initial_guess)
            initial_guess = TaskBasedProduction.find_initial_guess_gen(z, b_g, e_h, 'threshold', 1e-2, 'verbose', false);
        end
        [q, xT, fval] =  TaskBasedProduction.prod_fun_general(labor_input, z, b_g, e_h, 'initial_guess', initial_guess);
    end

    H = length(e_h);
    temp = zeros(H-1, 1); % Pre-allocate array for ratio values

    % Calculate the ratio e_{h} / e_{h-1} for h = 2:H and evaluate at xT[h-1]
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


