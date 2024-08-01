clear
clc
import TaskBasedProduction.*

% Initialize parameters
theta = 1.0;
kappa = 0.5;
z = 1.2;
alphaVec = [0.1, 0.2, 0.3];
labor_input=[0.5; 0.04; 0.19];

initial_guess=find_initial_guess(theta, kappa, z, alphaVec);
[q, xT, fval, initial_guess] = prod_fun(labor_input, theta, kappa, z, alphaVec, ...
    'initial_guess', initial_guess, ...
    'f_tol', 1e-4, ...  % Function tolerance for stopping criteria
    'iterations', 2000000); % Display output at each iteration

% Display the results
disp(['Quantity Produced: ', num2str(q)]);
disp(['Task Thresholds: ', num2str(xT')]); % Transpose xT for better display if it's a row vector
disp(['Approximation error: ', num2str(fval)]);

% Call unitInputDemand and print the output
labor_demand = unitInputDemand(xT, q, theta, kappa, z, alphaVec);
disp('Labor Demand:');
disp(labor_demand);


% Call margProdLabor with labor demand
mpl = margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q);

disp('Marginal Products of Labor:');
disp(mpl);



% Call elasticity_sub_comp
[epsilon_h_sub, epsilon_h_compl] = elasticity_sub_comp(labor_input, theta, kappa, z, alphaVec, mpl, xT, q);

disp('Elasticity of Substitution:');
disp(epsilon_h_sub);
disp('Elasticity of Complementarity:');
disp(epsilon_h_compl);

%% General parameterization

% Define b_g as a function handle
b_g = @(x) (x.^(kappa-1) .* exp(-x/theta)) / (theta^kappa * gamma(kappa));

% Define e_h functions as function handles
e_h1 = @(x) exp(0.1 * x);
e_h2 = @(x) exp(0.2 * x);
e_h3 = @(x) exp(0.3 * x);
e_h = {e_h1, e_h2, e_h3}; % Cell array of function handles
initial_guess = find_initial_guess_gen(z, b_g, e_h, 'threshold', 1e-2, 'verbose', false);
[q_gen, xT_gen, fval, initial_guess] =  prod_fun_general(labor_input, z, b_g, e_h, 'initial_guess', initial_guess);
% Display the results
disp(['Quantity Produced: ', num2str(q_gen)]);
disp(['Task Thresholds: ', num2str(xT_gen')]); % Transpose xT for better display if it's a row vector
disp(['Approximation error: ', num2str(fval)]);

% Call the unitInputDemand_general function
labor_demand_general = unitInputDemand_general(xT_gen, q_gen, z, b_g, e_h);

% Display the result
disp('Labor Demand:');
disp(labor_demand_general);

mpl_gen=margProdLabor_general(labor_demand_general, z, b_g, e_h, [],[], initial_guess);
[epsilon_h_sub_gen, epsilon_h_compl_gen] =  elasticity_sub_comp_general(labor_input, z, b_g, e_h, mpl_gen,[], []);
disp('Elasticity of Substitution:');
disp(epsilon_h_sub_gen);
disp('Elasticity of Complementarity:');
disp(epsilon_h_compl_gen);

