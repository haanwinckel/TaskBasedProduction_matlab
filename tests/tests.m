clearxbar
clc
% Add the package folder to the MATLAB path
addpath(genpath('C:\Users\lenovo\.julia\dev\TaskBasedProduction_matlab\TaskBasedProduction'));
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

initial_guess=find_initial_guess(theta, kappa, z, alphaVec);
[q, xT, fval, initial_guess] = prod_fun(labor_input, theta, kappa, z, alphaVec, 'initial_guess', initial_guess);

mpl = margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q);
[epsilon_h_sub, epsilon_h_compl] = elasticity_sub_comp(labor_input, theta, kappa, z, alphaVec, mpl, xT);

% Test for unitInputDemand
labor_input2 = q*unitInputDemand(xT,theta, kappa, z, alphaVec);
assert(isapprox(labor_input, labor_input2, 1e-5), 'unitInputDemand test failed');

% Test for MPL function
mpl2 = margProdLabor(labor_input, theta, kappa, z, alphaVec, [],[], initial_guess);
assert(isapprox(mpl, mpl2, 1e-5), 'MPL function test failed');

% Numerical comparison for MPL and elasticity of complementarity
H = length(labor_input);
tol = 1e-4;
num_mpl = zeros(H, 1);
epsilon_compl_numerical = zeros(H, H);
for h = 1:H
    perturbation = zeros(H, 1);
    perturbation(h) = tol;
    [qp, xTp, fvalp, initial_guess_p] = prod_fun(labor_input + perturbation, theta, kappa, z, alphaVec, 'initial_guess', initial_guess);
    [qn, xTn, fvaln, initial_guess_n] = prod_fun(labor_input - perturbation, theta, kappa, z, alphaVec, 'initial_guess', initial_guess);
    num_mpl(h) = (qp - qn) / (2 * tol);
    for hprime = h+1:H
            MPL_d = numerical_second_deriv(labor_input, theta, kappa, z, alphaVec, h, hprime, tol, xTp, xTn, qp, qn);
            epsilon_compl_numerical(h, hprime) = q * MPL_d / (mpl(h) * mpl(hprime));
    end
     
end
assert(isapprox(mpl, num_mpl, 1e-2), 'MPL numerical comparison test failed');

% Find initial guess for general case
initial_guess_gen = find_initial_guess_gen(z, b_g, e_h, 'threshold', 1e-2, 'verbose', false);
[q_gen, xT_gen, fval, initial_guess_gen] = prod_fun_general(labor_input, z, b_g, e_h, 'initial_guess', initial_guess_gen);
labor_input_general = q_gen * unitInputDemand_general(xT_gen, z, b_g, e_h);
mpl_gen = margProdLabor_general(labor_input_general, z, b_g, e_h, xT_gen, q_gen);
[epsilon_sub_gen, epsilon_compl_gen] = elasticity_sub_comp_general(labor_input_general, z, b_g, e_h, mpl_gen, xT_gen);

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
    unitInputDemand_general(xT_gen, z, b_g_invalid, e_h);
    error('Invalid density test failed');
catch
    disp('Invalid density test passed');
end

disp('ALL TESTS PASSED!!');

