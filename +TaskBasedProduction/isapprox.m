function result = isapprox(a, b, atol)
    % isapprox Checks if two arrays are approximately equal within a tolerance
    %
    % Usage:
    %   result = isapprox(a, b) % Uses default atol = 1e-5
    %   result = isapprox(a, b, atol) % Uses specified atol
    %
    % Inputs:
    %   a - First array
    %   b - Second array
    %   atol - (optional) Absolute tolerance, default is 1e-5
    %
    % Output:
    %   result - Boolean indicating if arrays are approximately equal

    if nargin < 3
        atol = 1e-5; % Default tolerance
    end

    result = all(abs(a - b) < atol);
end