clear
clc
import TaskBasedProduction.*

% Initialize parameters
theta = 1.0;
kappa = 0.5;
z = 1.2;
alphaVec = [0.1, 0.2, 0.3];
labor_input=[0.5; 0.04; 0.19];

% Call the prodFun function
[q, xT] = prodFun(labor_input, theta, kappa, z, alphaVec);
% Display the results
disp(['Quantity Produced: ', num2str(q)]);
disp(['Task Thresholds: ', num2str(xT')]); % Transpose xT for better display if it's a row vector

% Call unitInputDemand and print the output
labor_demand = unitInputDemand(xT, q, theta, kappa, z, alphaVec);
disp('Labor Demand:');
disp(labor_demand);


% Call margProdLabor with labor demand
mpl = margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q);

disp('Marginal Products of Labor:');
disp(mpl);



% Call elasticity_sub_comp
[epsilon_h_sub, epsilon_h_compl] = elasticitySubComp(labor_input, theta, kappa, z, alphaVec, mpl, xT, q);

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
initial_guess = getStartGuessGen_xT(z, b_g, e_h);
[q_gen, xT_gen] =  prodFunGeneral(labor_input, z, b_g, e_h);

% Display the results
disp(['Quantity Produced: ', num2str(q_gen)]);
disp(['Task Thresholds: ', num2str(xT_gen')]); % Transpose xT for better display if it's a row vector

% Call the unitInputDemand_general function
labor_demand_general = unitInputDemandGeneral(xT_gen, q_gen, z, b_g, e_h);

% Display the result
disp('Labor Demand:');
disp(labor_demand_general);

mpl_gen=margProdLaborGeneral(labor_demand_general, z, b_g, e_h, [],[]);
mpl_gen=margProdLaborGeneral([], z, b_g, e_h, xT_gen, q_gen);
[epsilon_h_sub_gen, epsilon_h_compl_gen] =  elasticitySubCompGeneral(labor_input, z, b_g, e_h, mpl_gen);
disp('Elasticity of Substitution:');
disp(epsilon_h_sub_gen);
disp('Elasticity of Complementarity:');
disp(epsilon_h_compl_gen);

%% Examples by use cases
%% 
disp('Use Case 1: Competitive Market with Functional Forms')
q = 1;
wage = [0.1; 0.2; 0.7];

% Compute thresholds xT using the functional form
diff_alpha = diff(alphaVec);
log_wage_ratios = log(wage(2:end) ./ wage(1:end-1));
xT = (1 ./ diff_alpha) .* log_wage_ratios;

% Calculate labor unit input requirements
labor_input_1 = unitInputDemand(xT, q, theta, kappa, z, alphaVec);
disp("Labor Input:");
disp(labor_input_1);


disp('Use Case 2: Given Labor Input, Use the Production Function to Obtain Total Production and Task Thresholds')
labor_input = [0.5; 0.04; 0.19];
[q, xT] = prodFun(labor_input, theta, kappa, z, alphaVec);

% General parameterization
[q, xT] = prodFunGeneral(labor_input, z, b_g, e_h);

disp('Use Case 3: Calculate the Elasticity of Complementarity and Substitution')
% Call elasticity_substitution with labor demand given
[epsilon_sub, epsilon_compl] = elasticitySubComp(labor_input, theta, kappa, z, alphaVec, [], [],[]);
disp("Allen partial elasticity of substitution:");
disp(epsilon_sub);
disp("Hicks partial elasticity of substitution:");
disp(epsilon_compl);
[epsilon_sub, epsilon_compl] = elasticitySubComp([], theta, kappa, z, alphaVec, mpl, xT, q)

disp('Use Case 4: Problem of the Firm in a Monopsonistic Labor Market')
% Define parameters
beta = 4;
L = [1; 1; 1];
p = 1;
w_inclusive = [0.4; 0.9; 2];
initial_guess = getStartGuess_xT(theta, kappa, z, alphaVec);

% Set optimization options for fmincon
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

% Perform the optimization using fmincon
% Pass the additional parameters using an anonymous function
result = fmincon(@(x) objective_to_minimize(x, theta, kappa, z, alphaVec, beta, L, p, w_inclusive), initial_guess, [], [], [], [], [], [], [], options);

% Extract optimized results
optimal_initial_guess = result;
q_opt = exp(optimal_initial_guess(1));
xT_opt = cumsum(exp(optimal_initial_guess(2:end)));

% Display results
disp('Optimal q:');
disp(q_opt);
disp('Optimal xT:');
disp(xT_opt);

%% General blueprint and efficiency functions
disp('Use case 1: Competitive labor market with general functions')
% Objective function for optimization
objective = @(x, h) (e_h{h+1}(x) / e_h{h}(x) - wage(h+1) / wage(h))^2;

% Find the solutions for xT_h
H = length(wage);
xT = NaN(H-1, 1);

% Initial guess: Call getStartGuessGen_xT to get sensible starting point
initial_guess = getStartGuessGen_xT(z, b_g, e_h);

for h = 1:H-1
    x0 = initial_guess(h+1);
    % Perform the optimization
    res = fminunc(@(x) objective(x, h), x0);
    % Extract the optimized value
    xT(h) = res;
end
disp('Use Case 2: Given Labor Input, Use the Production Function to Obtain Total Production and Task Thresholds')
% General parameterization
[q, xT] = prodFunGeneral(labor_input, z, b_g, e_h);


labor_input_2 = unitInputDemandGeneral(xT, q, z, b_g, e_h);
disp("Labor Input with general functions:");
disp(labor_input_2);
disp('Use Case 3: Calculate the Elasticity of Complementarity and Substitution')
[epsilon_h_sub_gen, epsilon_h_compl_gen] =  elasticitySubCompGeneral(labor_input, z, b_g, e_h, mpl_gen);
disp("Allen partial elasticity of substitution:");
disp(epsilon_h_sub_gen);
disp("Hicks partial elasticity of substitution:");
disp(epsilon_h_compl_gen);
disp('Use Case 4: Problem of the Firm in a Monopsonistic Labor Market (general blueprints and efficiency functions)')
% Perform the optimization using fmincon
result = fmincon(@(x) objective_to_minimize_general(x, z, b_g, e_h, beta, L, p, w_inclusive), ...
                 initial_guess, [], [], [], [], [], [], [], options);

% Extract optimized results
optimal_initial_guess = result;
q_opt = exp(optimal_initial_guess(1));
xT_opt = cumsum(exp(optimal_initial_guess(2:end)));

% Display results
disp('Optimal q:');
disp(q_opt);
disp('Optimal xT:');
disp(xT_opt);