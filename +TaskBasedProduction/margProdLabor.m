function mpl = margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q, initial_guess)
    % margProdLabor Calculates marginal products of labor
    %
    % Calculate marginal products of labor for each worker type given an array
    % of H labor demand values, an array of comparative advantage values alphaVec,
    % and optional arrays for task thresholds xT and quantity produced q. If xT or q
    % are not provided, they will be computed within the function.
    %
    % Inputs:
    %   labor_input - array of H labor demand values (vector)
    %   theta - blueprint scale parameter (scalar)
    %   kappa - blueprint shape parameter (scalar)
    %   z - productivity parameter (scalar)
    %   alphaVec - array of comparative advantage values (vector)
    %   xT - (optional) array of H-1 task thresholds (vector)
    %   q - (optional) scalar representing the quantity produced
    %   initial_guess - (optional) Initial guess for the task thresholds, used when computing xT and q.
    %
    % Output:
    %   mpl - marginal products of labor (vector)
    
    if nargin < 6 || isempty(xT) || isempty(q)
        if isempty(initial_guess)
            initial_guess = TaskBasedProduction.find_initial_guess(theta, kappa, z, alphaVec);  % Provide a default initial guess if not specified
        end
        [q, xT, fval] = TaskBasedProduction.prod_fun(labor_input, theta, kappa, z, alphaVec, 'initial_guess', initial_guess);
    end
    
    mpl_over_mpl1 = [1; cumprod(exp(diff(alphaVec') .* xT))];
    mpl1 = q / sum(mpl_over_mpl1 .* labor_input);
    mpl = mpl_over_mpl1 * mpl1;
end


