clear
clc
% Add the package folder to the MATLAB path
addpath(genpath('C:\Users\lenovo\.julia\dev\TaskBasedProduction_matlab'));
import TaskBasedProduction.*

% Initialize common parameters for tests
theta = 1.0;
kappa = 0.5;
z = 1.2;
alphaVec = [0.1, 0.2, 0.3];
labor_input=[0.5; 0.04; 0.19];

% Define the density function b_g(x)
b_g = @(x) (x.^(kappa-1) .* exp(-x/theta)) / (theta^kappa * gamma(kappa));
e_h1 = @(x) exp(0.1 * x);
e_h2 = @(x) exp(0.2 * x);
e_h3 = @(x) exp(0.3 * x);
e_h = {e_h1, e_h2, e_h3}; % Example e_h functions


[q, xT] = prodFun(labor_input, theta, kappa, z, alphaVec);

mpl = margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q);
[epsilon_h_sub, epsilon_h_compl] = elasticitySubComp(labor_input, theta, kappa, z, alphaVec, mpl, xT, q);

% Test for unitInputDemand
labor_input2 = unitInputDemand(xT, q, theta, kappa, z, alphaVec);
assert(isapprox(labor_input, labor_input2, 1e-5), 'unitInputDemand test failed');

% Test for MPL function
mpl2 = margProdLabor(labor_input, theta, kappa, z, alphaVec,xT, q);
assert(isapprox(mpl, mpl2, 1e-5), 'MPL function test failed');

% Numerical comparison for MPL and elasticity of complementarity
H = length(labor_input);
tol = 1e-4;
num_mpl = zeros(H, 1);
epsilon_compl_numerical = zeros(H, H);
for h = 1:H
    perturbation = zeros(H, 1);
    perturbation(h) = tol;
    [qp, xTp] = prodFun(labor_input + perturbation, theta, kappa, z, alphaVec);
    [qn, xTn] = prodFun(labor_input - perturbation, theta, kappa, z, alphaVec);
    num_mpl(h) = (qp - qn) / (2 * tol);
    for hprime = h+1:H
            MPL_d = numerical_second_deriv(labor_input, theta, kappa, z, alphaVec, h, hprime, xTp, xTn, qp, qn, []);
            epsilon_compl_numerical(h, hprime) = q * MPL_d / (mpl(h) * mpl(hprime));
    end
     
end
assert(isapprox(mpl, num_mpl, 1e-2), 'MPL numerical comparison test failed');

% Find initial guess for general case
[q_gen, xT_gen] = prodFunGeneral(labor_input, z, b_g, e_h);
labor_input_general = unitInputDemandGeneral(xT_gen, q_gen, z, b_g, e_h);
mpl_gen = margProdLaborGeneral(labor_input_general, z, b_g, e_h, xT_gen, q_gen);
[epsilon_sub_gen, epsilon_compl_gen] = elasticitySubCompGeneral(labor_input_general, z, b_g, e_h, mpl_gen, xT_gen);

% Test for prod_fun_general comparison
assert(isapprox(q_gen, q, 1e-5), 'prod_fun_general comparison test failed (q)');
assert(isapprox(xT_gen, xT, 1e-5), 'prod_fun_general comparison test failed (xT)');

% Test for labor input general
assert(isapprox(labor_input_general, labor_input, 1e-5), 'Labor input general test failed');

% Test for MPL general
assert(isapprox(mpl_gen, mpl, 1e-5), 'MPL general test failed');

% Test for elasticity of complementarity general
assert(all(isapprox(epsilon_h_compl, epsilon_compl_gen, 1e-2)), 'Elasticity compl general test failed');

% Test for elasticity of substitution general
assert(all(isapprox(epsilon_h_sub, epsilon_sub_gen, 20)), 'Elasticity sub general test failed');

% Test for invalid density
b_g_invalid = @(x) 4 * exp(-x); % This does not integrate to 1 over the entire domain
try
    unitInputDemandGeneral(xT_gen, z, b_g_invalid, e_h);
    error('Invalid density test failed');
catch
    disp('Invalid density test passed');
end

disp('ALL TESTS PASSED!!');

