function df = numerical_derivative(f, x, delta)
    if nargin < 3
        delta = 1e-5; % Default value for delta
    end
    df = (f(x + delta) - f(x - delta)) / (2 * delta); % Central difference method
end
