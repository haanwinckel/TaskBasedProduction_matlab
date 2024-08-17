function obj_value = objective_to_minimize(initial_guess, theta, kappa, z, alphaVec, beta, L, p, w_inclusive)
    % Compute q and xT from the initial guess
    q = exp(initial_guess(1));
    xT = cumsum(exp(initial_guess(2:end)));

    % Calculate labor input demand
    labor_input = TaskBasedProduction.unitInputDemand(xT, q, theta, kappa, z, alphaVec);
    
    % Calculate marginal product of labor (MPL)
    MPL = TaskBasedProduction.margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q);
    
    % Calculate wages
    w = p * (beta / (beta + 1)) * MPL;
    
    % Calculate labor supply
    labor_supply = (w ./ w_inclusive) .^ beta .* L;
    
    % Objective to minimize: sum of squared log differences between labor input and labor supply
    obj_value = sum(log(labor_input ./ labor_supply).^2);
end