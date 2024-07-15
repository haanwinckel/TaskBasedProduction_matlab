function epsilon_compl_numerical = numerical_elasticity_compl(labor_demand_specific, q, MPL,  theta, kappa, z, alphaVec, xT)
    H = length(MPL);
    epsilon_compl_numerical = zeros(H, H);
    
    for h = 1:H
        for hprime = h+1:H
            MPL_d = TaskBasedProduction.numerical_second_deriv(labor_demand_specific, theta, kappa, z, alphaVec, xT, h, hprime);
            epsilon_compl_numerical(h, hprime) = q * MPL_d / (MPL(h) * MPL(hprime));
        end
    end
end