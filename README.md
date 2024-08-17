# TaskBasedProduction
TaskBasedProduction is a MATLAB toolbox that provides functions for calculating unit labor demands, marginal products of labor, production functions, assignment thresholds, elasticities of substitution, and complementarities among worker types in a task-based production model. The package includes utilities for handling incomplete gamma functions and power series representations to facilitate these calculations. This package was developed by Daniel Haanwinckel and Luca Lorenzini based on the paper *Supply, Demand, Institutions, and Firms: A Theory of Labor Market Sorting and the Wage Distribution* by Daniel Haanwinckel ([NBER Working Paper No. 31318](https://www.nber.org/papers/w31318)).


## Installation
To install TaskBasedProduction, you can clone the repository and run the install.m script to add it to your MATLAB path:
```matlab
% Clone the repository
git clone https://github.com/yourusername/TaskBasedProduction.git

% Change directory to the cloned repository
cd TaskBasedProduction

% Run the install script
install

## Usage Example
```matlab
% Clear workspace and import package
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

mpl_gen=margProdLaborGeneral(labor_demand_general, z, b_g, e_h, [],[], initial_guess);
mpl_gen=margProdLaborGeneral([], z, b_g, e_h, xT_gen, q_gen);
[epsilon_h_sub_gen, epsilon_h_compl_gen] =  elasticitySubCompGeneral(labor_input, z, b_g, e_h, mpl_gen);
disp('Elasticity of Substitution:');
disp(epsilon_h_sub_gen);
disp('Elasticity of Complementarity:');
disp(epsilon_h_compl_gen);
```


## Use Cases
This package provides models and functions for analyzing labor markets, task thresholds, and elasticity in both competitive and monopsonistic settings. It includes functions for both the Exponential-gamma parameterization from the paper (Section 4.4 of the paper) and for arbitrary efficiency functions and blueprints. The first examples demonstrate the Exponential-gamma parameterization, with general examples for arbitrary efficiency functions and blueprints provided later.
## Exponential-gamma parameterization
### Parameters

```matlab
theta = 1.0;
kappa = 0.5;
z = 1.2;
alphaVec = [0.1, 0.2, 0.3];
```


## 1) Competitive labor market
In a competitive labor market, wages are taken as given and must equate the marginal product of labor. Optimality requires that marginal product ratios equal wage ratios, a relation that is used to obtain the task thresholds xT. Once the task thresholds are known, the function unitInputDemand is used to obtain the labor per output, setting q = 1. 
# Competitive labor market with functional forms
```matlab
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
```

## 2) Given labor input, use the production function to obtain total production and task thresholds
If labor inputs per each type are known, they can be given to the function prodFun or prodFunGeneral to compute the task thresholds and the total output produced.
### With functional forms
```matlab
labor_input = [0.5; 0.04; 0.19];
[q, xT] = prodFun(labor_input, theta, kappa, z, alphaVec);
```


## 3)  Elasticity of complementarity and substitution
Use the function elasticitySubComp to obtain the elasticities of substitution and complementarity. Precompiled values for marginal products (MPL), task thresholds (xT), and total output (q) can be used for efficiency, but they are computed within the function if not provided.

### With labor demand given
```matlab
[epsilon_sub, epsilon_compl] = elasticitySubComp(labor_input, theta, kappa, z, alphaVec, [], [],[]);
disp("Allen partial elasticity of substitution:");
disp(epsilon_sub);
disp("Hicks partial elasticity of substitution:");
disp(epsilon_compl);
```
### With Task Thresholds and Total Output
```matlab 
[epsilon_sub, epsilon_compl] = elasticitySubComp([], theta, kappa, z, alphaVec, mpl, xT, q)
disp("Allen partial elasticity of substitution:");
disp(epsilon_sub);
disp("Hicks partial elasticity of substitution:");
disp(epsilon_compl);
```

## 4) Problem of the firm in a monopsonistic labor market:
First, define the parameter governing the elasticity of labor supply:
``` matlab
beta = 4;
```
Then define general equilibrium objects (total labor supply and inclusive value) and the price of the good sold by the firm:
``` matlab
# Define general equilibrium objects
L = [1; 1; 1];
p = 1;
w_inclusive = [0.4; 0.9; 2];
``` 
Define an objective function that will be minimized to find the solution in total output (q) and task thresholds (xT). This function takes a guess for q and xT as input, computes the labor input required by the guess, and the marginal product of labor. It then computes the wage as a constant markdown applied to the marginal product of labor and the implied labor supply given the firm wage. Finally, the function returns an error that is the discrepancy from the labor supply to the labor required implied by the guess xT and q.

```  matlab
% Define the objective function for optimization derived from the firm problem 
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
``` 
To obtain the result, first a good initial guess is obtained and then the above function is minimized.

``` matlab
% Initial guess
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
``` 



## Arbitrary efficiency functions and blueprints
### Parameters

```matlab
% Define b_g as a function handle
b_g = @(x) (x.^(kappa-1) .* exp(-x/theta)) / (theta^kappa * gamma(kappa));

% Define e_h functions as function handles
e_h1 = @(x) exp(0.1 * x);
e_h2 = @(x) exp(0.2 * x);
e_h3 = @(x) exp(0.3 * x);
e_h = {e_h1, e_h2, e_h3}; % Cell array of function handles 
```


## 1) Competitive labor market
In a competitive labor market, wages are taken as given and must equate the marginal product of labor. Optimality requires that marginal product ratios equal wage ratios, a relation that is used to obtain the task thresholds xT. Once the task thresholds are known, the function unitInputDemand is used to obtain the labor per output, setting q = 1. 
```matlab
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

labor_input_2 = unitInputDemandGeneral(xT, q, z, b_g, e_h);
disp("Labor Input with general functions:");
disp(labor_input_2);
```

## 2) Given labor input, use the production function to obtain total production and task thresholds
If labor inputs per each type are known, they can be given to the function prod_fun or prod_fun_general to compute the task thresholds and the total output produced.
### With general parameterization
``` matlab
[q, xT] = prodFunGeneral(labor_input, z, b_g, e_h);
```

## 3)  Elasticity of complementarity and substitution
Use the function elasticity_sub_comp to obtain the elasticities of substitution and complementarity. Precompiled values for marginal products (MPL), task thresholds (xT), and total output (q) can be used for efficiency, but they are computed within the function if not provided.
```matlab
[epsilon_h_sub_gen, epsilon_h_compl_gen] =  elasticitySubCompGeneral(labor_input, z, b_g, e_h, mpl_gen);
```
## 4) Problem of the firm in a monopsonistic labor market:
First, define the parameter governing the elasticity of labor supply:
``` matlab
beta=4
```
Then define general equilibrium objects (total labor supply and inclusive value) and the price of the good sold by the firm:
``` matlab
# Define general equilibrium objects
L=[1 ; 1 ; 1] # Total labor force 
p=1  # Price for the good
w_inclusive=[0.4; 0.9; 2]  # Inclusive value of wages 
``` 
Define an objective function that will be minimized to find the solution in total output (q) and task thresholds (xT). This function takes a guess for q and xT as input, computes the labor input required by the guess, and the marginal product of labor. It then computes the wage as a constant markdown applied to the marginal product of labor and the implied labor supply given the firm wage. Finally, the function returns an error that is the discrepancy from the labor supply to the labor required implied by the guess xT and q.
``` matlab
% Same thing but with the general functions
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
```




## Functions and Features
### 1) **unitInputDemand**

Calculates unit labor demands given blueprint scale `theta`, blueprint shape `kappa`, productivity `z`, an array of comparative advantage values `alphaVec` with H elements (one for each worker type), and an array `xT` of H-1 thresholds in task space. The function handles both positive and negative `upsilon` values, corresponding to different conditions on labor demand.

The function also includes an optional parameter `skipParamChecks` that allows the user to bypass parameter validation for performance reasons.

### Arguments
- `xT`: A vector of H-1 thresholds in task space. Thresholds should be non-decreasing, and the function adds `0` and `Inf` to the threshold array for boundary conditions.
- `q`: The total quantity produced (scalar).
- `theta`: The blueprint scale parameter (scalar). Must be greater than 0.
- `kappa`: The blueprint shape parameter (scalar). Must be greater than 0.
- `z`: The productivity parameter (scalar). Must be greater than 0.
- `alphaVec`: A vector of comparative advantage values with H elements (one for each worker type). The values in `alphaVec` must be strictly increasing.
- `skipParamChecks` (optional): A boolean flag indicating whether to skip parameter checks (default is `false`). If `false`, the function validates the inputs to ensure they meet the expected criteria.

### Returns
- `labor_input`: A vector representing the labor demand for each worker type. The length of this vector is H, where H is the number of labor types (derived from the length of `xT` plus one).

The labor demands are computed based on the values of `upsilon`, which depends on `alphaVec` and `theta`. The function handles different cases for positive, zero, and negative `upsilon` values, using helper functions to calculate the required labor demand under each condition.

### Example Usage
```matlab
labor_input = unitInputDemand(xT, q, theta, kappa, z, alphaVec, skipParamChecks);
 ```

### 2) **getStartGuess_xT**

Generates an initial guess for the optimization problem in `prod_fun` such that the implied labor demand is non-trivial.

This function creates a starting guess for the production quantity `q` and task thresholds `xT` by sampling from a gamma distribution. The thresholds are adjusted iteratively to ensure that the implied labor demand for each task is above a specified threshold.

### Arguments
- `theta`: The scale parameter of the gamma distribution.
- `kappa`: The shape parameter of the gamma distribution.
- `z`: A scaling factor for the labor input.
- `alphaVec`: An array of task-specific parameters with H elements (one for each task).
- `threshold` (optional): The minimum acceptable labor demand for each task (default is `1e-2`).

### Returns
- `initial_guess`: A vector containing the initial guess for the optimization, which includes the log of the initial production quantity `q` (starting with `log(1) = 0`) and the initial task thresholds `xT`.

### Example Usage
```matlab
initial_guess = getStartGuess_xT(theta, kappa, z, alphaVec, threshold);
```

 ### 3) **prodFun**
Calculates the quantity produced (`q`), and task thresholds (`xT`) given labor inputs (`labor_input`), blueprint scale `theta`, blueprint shape `kappa`, productivity `z`, and an array of comparative advantage values `alphaVec` with H elements (one for each worker type).

**Inputs:**
- `labor_input`: Array of labor inputs of different types.
- `theta`: Blueprint scale parameter.
- `kappa`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `alphaVec`: Array of comparative advantage values with H elements.
- `initial_guess`: (optional) Initial guess for optimization. If not provided, it is generated using `getStartGuess_xT`.
- `x_tol`: (optional) Tolerance for x values. Default is the optimization's default.
- `f_tol`: (optional) Tolerance for function values. Default is the optimization's default.
- `g_tol`: (optional) Tolerance for gradient values. Default is the optimization's default.
- `iterations`: (optional) Maximum number of iterations. Default is 1000000.
- `max_retries`: (optional) Maximum number of retries. Default is 5.
- `verbose`: (optional) If `true`, displays detailed information during optimization. Default is `false`.

**Returns:**
- `q`: Quantity produced (scalar).
- `xT`: Array of task thresholds (vector).
- `fval`: Final value of the objective function.
- `initial_guess`: Initial guess for the optimization that worked well.
``` matlab
     prodFun(labor_input, theta, kappa, z, alphaVec, varargin)
```
### 4) **margProdLabor**
Calculates the marginal productivity of labor for each worker type given the input parameters.

**Arguments:**
- `labor_input`: An array of labor demand values.
- `theta`: Blueprint scale parameter.
- `kappa`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `alphaVec`: An array of comparative advantage values.
- `xT`: (optional) An array representing the precomputed task thresholds. If not provided, it will be computed within the function.
- `q`: (optional) A scalar representing the precomputed quantity produced. If not provided, it will be computed within the function.
- `initial_guess`: (optional) A vector containing the precomputed initial guess for the optimization to compute xT if not provided. If not provided, it will be computed within the function.

**Returns:**
- An array representing the marginal productivity of labor for each worker type.

```matlab
    margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q)
```

### 5) **elasticitySubComp**

Calculates the elasticity of substitution and complementarity for a given set of parameters.

**Arguments:**
- `labor_input`: An array of labor inputs of different types with H elements.
- `theta`: Blueprint scale parameter.
- `kappa`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `alphaVec`: An array of comparative advantage values with H elements.
- `MPL`: (optional) An array representing the marginal productivity of labor. If not provided, it will be computed within the function.
- `xT`: (optional) An array representing precomputed task thresholds. If not provided, it will be computed within the function.
- `initial_guess`: (optional) A vector containing the precomputed initial guess for the optimization to compute xT if not provided. If not provided, it will be computed within the function.

**Returns:**
- `ϵ_h_sub`: Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns).
- `ϵ_h_compl`: Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns).

``` matlab
    elasticitySubComp(labor_input, theta, kappa, z, alphaVec, MPL, xT, q)
```
### 6) **unitInputDemandGeneral**

Calculates unit labor demands given an array `xT` of H-1 thresholds in task space, a productivity value `z`, 
a density function `b_g` for the task distribution, and an array `e_h` of H functions representing the cost of each labor type as a function of task complexity.

**Arguments:**
- `xT`: A vector of H-1 thresholds in task space.
- `z`: Productivity value.
- `b_g`: A density function for the task distribution.
- `e_h`: A vector of H functions representing the cost of each labor type as a function of task complexity.
- `q`: Total production.

**Returns:**
- A vector representing the labor demand for each labor type.

``` matlab 
unitInputDemandGeneral(xT, q, z, b_g, e_h)
```

### 7) **getStartGuessGen_xT**

Generates an initial guess for the optimization problem using a general density function such that the implied labor demand is non-trivial.

**Arguments:**
- `z`: A scaling factor for the labor input.
- `b_g`: The general density function.
- `e_h`: A vector of task-specific functions.
- `threshold`: (optional) The minimum acceptable labor demand for each task (default = 1e-2).
- `verbose`: (optional) If true, display detailed output information (default = false).

**Returns:**
- `initial_guess`: A vector containing the initial guess for the optimization, including the log of the initial production quantity `q` and the initial task thresholds `xT`.

``` matlab
getStartGuessGen_xT(z, b_g, e_h, varargin)
``` 


# 8) **prodFunGeneral**

Calculates the quantity produced (q), and task thresholds (xT) given labor inputs (labor_input), productivity z, a general blueprint density function (b_g), and a vector of efficiency functions (e_h) for each labor type.

## Inputs:
- `labor_input`: Array of labor inputs of different types.
- `z`: Productivity parameter.
- `b_g`: Blueprint density function.
- `e_h`: Vector of efficiency functions, one for each type.
- `initial_guess`: (optional) Initial guess for optimization. If not provided, it will be generated.
- `x_tol`: (optional) Tolerance for the solution vector. Default is 1e-6.
- `f_tol`: (optional) Tolerance for the function value. Default is 1e-6.
- `g_tol`: (optional) Tolerance for the gradient. Default is 1e-6.
- `iterations`: (optional) Maximum number of iterations for the optimization. Default is 1000.
- `max_retries`: (optional) Maximum number of retries if the optimization fails. Default is 50.
- `verbose`: (optional) If true, display detailed information during optimization. Default is false.

## Returns:
- `q`: Quantity produced.
- `xT`: Array of task thresholds.

``` matlab
  prodFunGeneral(labor_input, z, b_g, e_h, varargin)
```


# 9) **margProdLaborGeneral**

Calculates the marginal productivity of labor for each worker type given the input parameters.

## Arguments
- `labor_input`: An array of labor inputs of different types with H elements.
- `z`: A productivity scalar.
- `b_g`: A task density function.
- `e_h`: A vector of comparative advantage functions.
- `xT`: (optional) An array representing the precomputed task thresholds. If not provided, it will be computed within the function.
- `q`: (optional) A scalar representing the precomputed quantity produced. If not provided, it will be computed within the function.
- `initial_guess`: (optional) A vector containing the precomputed initial guess for the optimization to compute xT if not provided. If not provided, it will be computed within the function.

## Returns
- `mpl`: An array representing the marginal productivity of labor for each worker type.

```matlab
function mpl = margProdLaborGeneral(labor_input, z, b_g, e_h, xT, q, initial_guess)
```

# 10) **elasticitySubCompGeneral**

Calculates the elasticity of substitution and complementarity for a given set of parameters.

## Arguments
- `labor_input`: An array of labor inputs of different types with H elements.
- `z`: Productivity parameter.
- `b_g`: General task density function.
- `e_h`: Vector of comparative advantage functions.
- `MPL`: (optional) An array representing the marginal productivity of labor. If not provided, it will be computed within the function.
- `xT`: (optional) An array representing precomputed task thresholds. If not provided, it will be computed within the function.
- `q`: (optional) A scalar representing total production. If not provided, it will be computed within the function.
- `initial_guess`: (optional) A vector containing the precomputed initial guess for the optimization to compute xT if not provided. If not provided, it will be computed within the function.

## Returns
- `ϵ_h_sub`: Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns).
- `ϵ_h_compl`: Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns).

```matlab
function [epsilon_h_sub, epsilon_h_compl] = elasticitySubCompGeneral(labor_input, z, b_g, e_h, MPL, xT, q)
```

## Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue if you have any suggestions or find any bugs.
