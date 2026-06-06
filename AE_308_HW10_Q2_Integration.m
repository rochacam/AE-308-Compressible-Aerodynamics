% Code Name: Capsule, numerical integration
% Code Description: This code will numerical integrate the afound
% analytical solution for a capsule like section for Cd in hypersonic flow.
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE308 - Section 1
% Date: 11/12/2025
% Worked With: N/a

%f(x) = cos^2(-x / sqrt(9 - x^2))
% over the interval from x = 0 to x = -2.

% 1. Define the limits of integration
a = 0;   % Lower limit
b = -2;  % Upper limit

% 2. Define the integrand function handle
% Note: The argument of the cosine function is treated as radians.
% The use of '^2' for squaring is vectorized with 'cos(x).^2'
integrand = @(x) cos( (-x) ./ sqrt(9 - x.^2) ).^2;

% 3. Perform the numerical integration
% The 'integral' function handles the integration from the lower limit (a) 
% to the upper limit (b), even when a > b.
  I = integral(integrand, a, b);
    
    % 4. Display the result
    fprintf('The numerical value of the integral is:\n');
    fprintf('I = %.6f\n', I);