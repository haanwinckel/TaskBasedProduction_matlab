function obj_value = objective_to_minimize_general(initial_guess, z, b_g, e_h, beta, L, p, w_inclusive)
    % Compute q and xT from the initial guess
    q = exp(initial_guess(1));
    xT = cumsum(exp(initial_guess(2:end)));
    
    % Calculate labor input demand using general functions
    labor_input = TaskBasedProduction.unitInputDemandGeneral(xT, q, z, b_g, e_h);
    
    % Calculate marginal product of labor (MPL)
    MPL = TaskBasedProduction.margProdLaborGeneral(labor_input, z, b_g, e_h, xT, q);
    
    % Calculate wages
    w = p * (beta / (beta + 1)) * MPL;
    
    % Calculate labor supply
    labor_supply = (w ./ w_inclusive) .^ beta .* L;
    
    % Objective to minimize: sum of squared log differences between labor input and labor supply
    obj_value = sum(log(labor_input ./ labor_supply).^2);
end