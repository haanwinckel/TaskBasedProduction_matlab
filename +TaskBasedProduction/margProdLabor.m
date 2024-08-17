function mpl = margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q)
    % margProdLabor Calculates marginal products of labor for each worker type
    %
    % Arguments:
    %   labor_input - array of labor demand values (vector). If empty, it will be computed internally.
    %   theta - blueprint scale parameter (scalar)
    %   kappa - blueprint shape parameter (scalar)
    %   z - productivity parameter (scalar)
    %   alphaVec - array of comparative advantage values (vector)
    %   xT - (optional) array representing the task thresholds (vector)
    %   q - (optional) scalar representing the quantity produced
    %
    % Returns:
    %   mpl - marginal products of labor (vector)
    
    % Check if xT or q are missing; if so, compute them using prodFun
    if nargin < 6 || isempty(xT) || isempty(q)
        [q, xT] = TaskBasedProduction.prodFun(labor_input, theta, kappa, z, alphaVec);
    end
    
    % If labor_input is missing, compute it using unitInputDemand
    if isempty(labor_input)
        labor_input = TaskBasedProduction.unitInputDemand(xT, q, theta, kappa, z, alphaVec);
    end
    
    % Compute mpl_over_mpl1 using the cumulative product of the exponentiated differences
    mpl_over_mpl1 = [1; cumprod(exp(diff(alphaVec') .* xT))];
    
    % Calculate mpl1
    mpl1 = q / sum(mpl_over_mpl1 .* labor_input);
    
    % Compute final marginal product of labor (mpl)
    mpl = mpl_over_mpl1 * mpl1;
end



