function MPL_d = numerical_second_deriv(labor_demand_specific, theta, kappa, z, alphaVec, xT, h, hprime, hstep)
    if nargin < 10
        hstep = 1e-4;
    end
    perturbation = zeros(length(xT)+1, 1);
    perturbation(h) = hstep;
    [qp, xTp] = TaskBasedProduction.prod_fun(labor_demand_specific + perturbation, theta, kappa, z, alphaVec);
    [qpp, xTpp] = TaskBasedProduction.prod_fun(labor_demand_specific - perturbation, theta, kappa, z, alphaVec);
    
    MPL_plus_h = TaskBasedProduction.margProdLabor2(theta, kappa, z, alphaVec, xTp);
    MPL_minus_h = TaskBasedProduction.margProdLabor2(theta, kappa, z, alphaVec, xTpp);
    
    MPL_d = (MPL_plus_h(hprime) - MPL_minus_h(hprime)) / (2 * hstep);
end