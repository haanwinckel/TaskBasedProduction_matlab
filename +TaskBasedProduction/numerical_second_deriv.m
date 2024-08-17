% Numerical Elasticity of complementarity
function derivative = numerical_second_deriv(labor_input, theta, kappa, z, alphaVec, h, hprime, xTp, xTn, qp, qn,hstep)
    if isempty(hstep)
        hstep = 1e-4;
    end
    assert(h > 0, 'h must be a natural number (positive integer)');
    assert(hprime > 0, 'hprime must be a natural number (positive integer)');

    perturbation = zeros(length(labor_input), 1);
    perturbation(h) = hstep;

    MPL_plus_h = TaskBasedProduction.margProdLabor(labor_input + perturbation, theta, kappa, z, alphaVec, xTp, qp);
    MPL_minus_h = TaskBasedProduction.margProdLabor(labor_input - perturbation, theta, kappa, z, alphaVec, xTn, qn);

    derivative = (MPL_plus_h(hprime) - MPL_minus_h(hprime)) / (2 * hstep);
end