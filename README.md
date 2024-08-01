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
clear;
clc;
import TaskBasedProduction.*;

% Initialize parameters
theta = 1.0;
kappa = 0.5;
z = 1.2;
alphaVec = [0.1, 0.2, 0.3];
labor_input = [0.5; 0.04; 0.19];

% Find initial guess
initial_guess = find_initial_guess(theta, kappa, z, alphaVec);

% Calculate production function
[q, xT, fval, initial_guess] = prod_fun(labor_input, theta, kappa, z, alphaVec, 'initial_guess', initial_guess);

% Display the results
disp(['Quantity Produced: ', num2str(q)]);
disp(['Task Thresholds: ', num2str(xT')]); % Transpose xT for better display if it's a row vector
disp(['Approximation error: ', num2str(fval)]);

% Calculate unit labor demand
labor_demand = unitInputDemand(xT,q, theta, kappa, z, alphaVec);
disp('Labor Demand:');
disp(labor_demand);

% Calculate marginal products of labor
mpl = margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q);
disp('Marginal Products of Labor:');
disp(mpl);

% Calculate elasticity of substitution and complementarity
[epsilon_h_sub, epsilon_h_compl] = elasticity_sub_comp(labor_input, theta, kappa, z, alphaVec, mpl, xT, initial_guess);
disp('Elasticity of Substitution:');
disp(epsilon_h_sub);
disp('Elasticity of Complementarity:');
disp(epsilon_h_compl);

% General parameterization example
% Define b_g as a function handle
b_g = @(x) (x.^(kappa-1) .* exp(-x/theta)) / (theta^kappa * gamma(kappa));

% Define e_h functions as function handles
e_h1 = @(x) exp(0.1 * x);
e_h2 = @(x) exp(0.2 * x);
e_h3 = @(x) exp(0.3 * x);
e_h = {e_h1, e_h2, e_h3}; % Cell array of function handles

% Find initial guess for general parameterization
initial_guess = find_initial_guess_gen(z, b_g, e_h, 'threshold', 1e-2, 'verbose', false);

% Calculate production function with general parameterization
[q_gen, xT_gen, fval, initial_guess] = prod_fun_general(labor_input, z, b_g, e_h, 'initial_guess', initial_guess);

% Display the results
disp(['Quantity Produced: ', num2str(q_gen)]);
disp(['Task Thresholds: ', num2str(xT_gen')]); % Transpose xT for better display if it's a row vector
disp(['Approximation error: ', num2str(fval)]);

% Calculate unit labor demand for general parameterization
labor_demand_general =  unitInputDemand_general(xT_gen, q_gen, z, b_g, e_h);
disp('Labor Demand:');
disp(labor_demand_general);

% Calculate marginal products of labor for general parameterization
mpl_gen = margProdLabor_general(labor_demand_general, z, b_g, e_h, [], [], initial_guess);

% Calculate elasticity of substitution and complementarity for general parameterization
[epsilon_h_sub_gen, epsilon_h_compl_gen] = elasticity_sub_comp_general(labor_input, z, b_g, e_h, mpl_gen, xT_gen);
disp('Elasticity of Substitution:');
disp(epsilon_h_sub_gen);
disp('Elasticity of Complementarity:');
disp(epsilon_h_compl_gen);
```
## Functions and Features
# 1) **unitInputDemand**
Calculates unit labor demands given blueprint scale `theta`, blueprint shape `kappa`, productivity `z`, an array of comparative advantage values `alphaVec` with H elements (one for each worker type), and an array `xT` of H-1 thresholds in task space.

# Arguments
- `xT`: An array of H-1 thresholds in task space.
- `theta`: Blueprint scale parameter.
- `kappa`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `alphaVec`: An array of comparative advantage values with H elements.
- `skipParamChecks`: A boolean indicating whether to skip parameter checks (default is false).

# Returns
- An array representing the labor demand for each labor type.
```matlab
    unitInputDemand(xT,theta, kappa, z, alphaVec, skipParamChecks)
 ```

# 2) **find_initial_guess**

Generate an initial guess for the optimization problem in `prod_fun` such that the implied labor demand is non-trivial.

# Arguments
- `labor_input`: The observed labor input for each task.
- `theta`: The scale parameter of the gamma distribution.
- `kappa`: The shape parameter of the gamma distribution.
- `z`: A scaling factor for the labor input.
- `alphaVec`: An array of task-specific parameters.
- `threshold`: The minimum acceptable labor demand for each task.

# Returns
- `initial_guess`: A vector containing the initial guess for the optimization, including the log of the initial production quantity `q` and the initial task thresholds `xT`.



``` matlab
    find_initial_guess(theta, kappa, z, alphaVec, threshold)
```


 # 3) **prod_fun**
Calculates the quantity produced (q), and task thresholds (xT) given labor inputs (l), blueprint scale theta, blueprint shape kappa, productivity z, and an array of 
comparative advantage values alphaVec with H elements (one for each worker type).

Inputs:
- `labor_input`: Array of labor inputs of different types.
- `theta`: Blueprint scale parameter.
- `kappa`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `alphaVec`: Array of comparative advantage values with H elements.
- `initial_guess`: (optional) Initial guess for optimization. If not provided, defaults to zeros array.
- `optim_options`: (optional) Optimization options. If not provided, defaults to high tolerance values.
- `verbose`: (optional)  If true, display detailed information during optimization. Default is false.
Returns:
- `q`: Quantity produced.
- `xT`: Array of task thresholds.
- `fval`: Final value of the objective function.
- `initial_guess`: Initial guess for the optimization that worked well in the optimization .

``` matlab
     prod_fun(labor_input, theta, kappa, z, alphaVec, varargin)
```
# 4) **margProdLabor**
Calculates the marginal productivity of labor for each worker type given the input parameters.

# Arguments
- `labor_input`: An array of labor demand values.
- `theta`: Blueprint scale parameter.
- `kappa`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `alphaVec`: An array of comparative advantage values.
- `xT`: (optional) An array representing the precomputed task thresholds. If not provided, it will be computed within the function.
- `q`: (optional) A scalar representing the precomputed quantity produced. If not provided, it will be computed within the function.
- `initial_guess`: (optional) A vector containing the precomputed initial guess for the optimization to compute xT if not provided. If not provided, it will be computed within the function.
# Returns
- An array representing the marginal productivity of labor for each worker type.
```matlab
     margProdLabor(labor_input, theta, kappa, z, alphaVec, xT, q, initial_guess)
```

# 5) **elasticity_sub_comp**

Calculates the elasticity of substitution and complementarity for a given set of parameters.

# Arguments
- `labor_input`: An array of labor inputs of different types with H elements.
- `theta`: Blueprint scale parameter.
- `kappa`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `alphaVec`: An array of comparative advantage values with H elements.
- `MPL`: (optional) An array representing the marginal productivity of labor. If not provided, it will be computed within the function.
- `xT`: (optional) An array representing precomputed task thresholds. If not provided, it will be computed within the function.
- `initial_guess`: (optional) A vector containing the precomputed initial guess for the optimization to compute xT if not provided. If not provided, it will be computed within the function.
# Returns
- `ϵ_h_sub`: Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns).
- `ϵ_h_compl`: Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns).

``` matlab
    elasticity_sub_comp(labor_input, theta, kappa, z, alphaVec, MPL, xT, initial_guess)
```
# 6) **unitInputDemand_general**

Calculates unit labor demands given an array `xT` of H-1 thresholds in task space, a productivity value `z`, 
a density function `b_g` for the task distribution, and an array `e_h` of H functions
representing the cost of each labor type as a function of task complexity.

The function first verifies that `b_g` is a valid density function. Then it computes
the labor demand for each labor type by numerically integrating the ratio `b_g(x) / (z * e_h[h](x))`
over the intervals defined by the thresholds in `xT`.

# Arguments
- `xT`: A vector of H-1 thresholds in task space.
- `z`: Productivity value.
- `b_g`: A density function for the task distribution.
- `e_h`: A vector of H functions representing the cost of each labor type as a function of task complexity.

# Returns
- A vector representing the labor demand for each labor type.


``` matlab 
unitInputDemand_general(xT, z, b_g, e_h)
```

# 7) **find_initial_guess_gen**

  
Generate an initial guess for the optimization problem using a general density function such that the implied labor demand is non-trivial.

# Arguments
- `labor_input`: The observed labor input for each task.
- `z`: A scaling factor for the labor input.
- `alphaVec`: An array of task-specific parameters.
- `pdf`: The general density function.
- `threshold`: The minimum acceptable labor demand for each task.
- `verbose`: (optional)  If true, display detailed output information. Default is false.
# Returns
- `initial_guess`: A vector containing the initial guess for the optimization, including the log of the initial production quantity `q` and the initial task thresholds `xT`.

``` matlab
find_initial_guess_gen(z, b_g, e_h, varargin)
``` 


# 8) **prod_fun_general**

Calculates the quantity produced (q), and task thresholds (xT) given labor inputs (labor_input), productivity z, general blueprint density function (b_g), and a vector of efficiency functions (e_h), one for each labor type.

Inputs:
- `labor_input`: Array of labor inputs of different types.
- `z`: Productivity parameter.
- `b_g`: Blueprint density function.
- `e_h`: Vector of efficiency functions, one for each type.
- `initial_guess`: (optional) Initial guess for optimization. If not provided, defaults to zeros array.
- `x_tol`: (optional) Tolerance for the solution vector. Default is [].
- `f_tol`: (optional) Tolerance for the function value. Default is 1e-6.
- `g_tol`: (optional) Tolerance for the gradient. Default is []].
- `iterations`: (optional) Maximum number of iterations for the optimization. Default is 1000000.
- `max_retries`: (optional) Maximum number of retries if the optimization fails. Default is 1000.
- `verbose`: (optional)  If true, display detailed information during optimization. Default is false.
Returns:
- `q`: Quantity produced.
- `xT`: Array of task thresholds.
- `fval`: Final value of the objective function.

``` matlab
 prod_fun_general(labor_input, z, b_g, e_h, varargin)
```


# 9) **margProdLabor_general**

Calculates the marginal productivity of labor for each worker type given the input parameters.

# Arguments
- `labor_input`: An array of labor inputs of different types with H elements.
- `z`: A productivity scalar.
- `b_g`: A task density function.
- `e_h`: A vector of comparative advantage functions.
- `xT`: (optional) An array representing the precomputed task thresholds. If not provided, it will be computed within the function.
- `q`: (optional) A scalar representing the precomputed quantity produced. If not provided, it will be computed within the function.
- `initial_guess`: (optional) A vector containing the precomputed initial guess for the optimization to compute xT if not provided. If not provided, it will be computed within the function.
# Returns
- An array representing the marginal productivity of labor for each worker type.

``` matlab
margProdLabor_general(labor_input, z, b_g, e_h, xT, q, initial_guess)
```

# 10) **elasticity_sub_comp_general**

Calculates the elasticity of substitution and complementarity for a given set of parameters.

# Arguments
- `labor_input`: An array of labor inputs of different types with H elements.
- `z`: Productivity parameter.
- `b_g`: General task density function.
- `e_h`: Vector of comparative advantage functions.
- `MPL`: (optional) An array representing the marginal productivity of labor. If not provided, it will be computed within the function.
- `xT`: (optional) An array representing precomputed task thresholds. If not provided, it will be computed within the function.
- `initial_guess`: (optional) A vector containing the precomputed initial guess for the optimization to compute xT if not provided. If not provided, it will be computed within the function.
# Returns
- `ϵ_h_sub`: Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns).
- `ϵ_h_compl`: Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns).

``` matlab
elasticity_sub_comp_general(labor_input, z, b_g, e_h, MPL, xT, initial_guess)
```

## Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue if you have any suggestions or find any bugs.
