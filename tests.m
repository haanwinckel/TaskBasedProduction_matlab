clear
clc

%% Initialize parameters for the first test
theta = 1.0;
kappa = 0.5;
z = 1.2;
alphaVec = [0.1, 0.2, 0.3];
xT = [0.4, 0.5];

% Expected values
expected_labor_demand = [0.5178877620105745; 0.04097664559663147; 0.18579978864645658];

% Call unitInputDemand and assert the output
labor_demand = TaskBasedProduction.unitInputDemand(theta, kappa, z, alphaVec, xT);
disp('Labor Demand:');
disp(labor_demand);
assert(all(abs(labor_demand(:) - expected_labor_demand(:)) < 1e-3), 'unitInputDemand test failed');

%% Test for margProdLabor with labor demand
expected_marginal_products = [1.309184899610229, 1.3626137489243066, 1.4324764497687001];
marginal_products = TaskBasedProduction.margProdLabor(labor_demand, alphaVec, xT);
disp('Marginal Products of Labor:');
disp(marginal_products);
assert(all(abs(marginal_products(:) - expected_marginal_products(:)) < 1e-5), 'margProdLabor with labor demand test failed');

%% Test for margProdLabor with blueprint characteristics
expected_marginal_products_alt = [1.309184899610229, 1.3626137489243066, 1.4324764497687001];
marginal_products_alt = TaskBasedProduction.margProdLabor2(theta, kappa, z, alphaVec, xT);
disp('Marginal Products of Labor (using blueprint characteristics):');
disp(marginal_products_alt);
assert(all(abs(marginal_products_alt(:) - expected_marginal_products_alt(:)) < 1e-3), 'margProdLabor with blueprint characteristics test failed');

%% Define b_g and e_h functions
b_g = @(x) (x.^(kappa-1) .* exp(-x/theta)) / (theta^kappa * gamma(kappa));
e_h1 = @(x) exp(0.1 * x);
e_h2 = @(x) exp(0.2 * x);
e_h3 = @(x) exp(0.3 * x);
e_h = {e_h1, e_h2, e_h3}; % Cell array of function handles

%% Test for unitInputDemand_general and comparison with unitInputDemand
% Call prod_fun with labor demand
[q, xT] = TaskBasedProduction.prod_fun(labor_demand, theta, kappa, z, alphaVec);
labor_demand_general = TaskBasedProduction.unitInputDemand_general(xT, z, b_g, e_h);
disp('General Labor Demand:');
disp(labor_demand_general);

[q_gen, xT_gen, fval] = TaskBasedProduction.prod_fun_general(labor_demand_general, z, b_g, e_h);
mpl_gen = TaskBasedProduction.margProdLabor_general(xT_gen, labor_demand_general, e_h);


disp('Quantity Produced:');
disp(q);
disp('Task Thresholds:');
disp(xT);
mpl=marginal_products;
% Call elasticity_sub_comp
[epsilon_h_sub, epsilon_h_compl] = TaskBasedProduction.elasticity_sub_comp(xT, labor_demand, q, marginal_products, theta, kappa, z, alphaVec);

disp('Elasticity of Substitution:');
disp(epsilon_h_sub);
disp('Elasticity of Complementarity:');
disp(epsilon_h_compl);

% Call elasticity_sub_comp_general
[epsilon_h_sub_gen, epsilon_h_compl_gen] = TaskBasedProduction.elasticity_sub_comp_general(xT_gen, labor_demand_general, q_gen, mpl_gen, z, b_g, e_h);

disp('General Elasticity of Substitution:');
disp(epsilon_h_sub_gen);
disp('General Elasticity of Complementarity:');
disp(epsilon_h_compl_gen);

% Assertions for comparison tests
assert(all(abs(labor_demand_general(:) - labor_demand(:)) < 1e-2), 'unitInputDemand comparison test failed');
assert(abs(q(:) - q_gen(:)) < 1e-2, 'Quantity comparison test failed');
assert(all(abs(xT(:)- xT_gen(:)) < 1e-2), 'Task thresholds comparison test failed');
assert(all(abs(mpl(:) - mpl_gen(:)) < 1e-2), 'MPL comparison test failed');
assert(all(abs(epsilon_h_sub(:) - epsilon_h_sub_gen(:)) < 30), 'Elasticity of substitution comparison test failed');
assert(all(abs(epsilon_h_compl(:) - epsilon_h_compl_gen(:)) < 1e-2), 'Elasticity of complementarity comparison test failed');

% Test for invalid density function
b_g_invalid = @(x) 4 * exp(-x); % This does not integrate to 1 over the entire domain
try
    TaskBasedProduction.unitInputDemand_general(xT, z, b_g_invalid, e_h);
    error('unitInputDemand_general did not throw an error for invalid density function');
catch
    disp('unitInputDemand_general threw an error for invalid density function as expected');
end

%% Numerical second derivative function
epsilon_compl2 = TaskBasedProduction.numerical_elasticity_compl(labor_demand, q, mpl, theta, kappa, z, alphaVec, xT);
assert(all(abs(epsilon_h_compl(:) - epsilon_compl2(:)) < 1e-3), 'Elasticity of complementarity test (epsilon_compl) failed');
assert(all(abs(epsilon_h_compl_gen(:) - epsilon_compl2(:)) < 1e-3), 'Elasticity of complementarity test (epsilon_compl_gen) failed');

