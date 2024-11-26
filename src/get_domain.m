function domain = get_domain(a, L, H, c, acc)
% get_domain - function to automatically generate a complex domain suitable
% to perform the needed calculations and plots for the two-airfoils
% configuration. Upper and lower boundaries for the real and imaginary
% coordinates of the complex plane are computed by inverse trsnsforming the
% points at x_max = L + 3*c, x_min = -3*c, y_max = H + 3*c, y_min = -3*c in
% the Joukowski profile plane.
% 
% INPUTS: 
% - a, float: a parameter of the Joukowski circle
% - L, float: x-coodinate of the tail LE (ref system has wing LE at (0, 0))
% - H, float: y-coodinate of the tail LE (ref system has wing LE at (0, 0))
% - c, float: tail chord
% - acc, float: accuracy of the discretization, default value 1000
%
% OUTPUTS:
% - domain, complex float: complex domain for the Joukowski circle
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 26/11/24, Boscariol Jacopo
%               Changes: release.
    
    % checking accuracy specification
    if nargin < 5
        acc = 1000;
    end

    % getting upper and lower domain boundaries from inverse transform
    up_real = (L + 3*c + sqrt((L + 3*c)^2 - 4*a^2))/2;
    up_imag = (H + 3*c + sqrt((H + 3*c)^2 - 4*a^2))/2;

    low_real = (-3*c - sqrt((-3*c)^2 - 4*a^2))/2;
    low_imag = (-3*c - sqrt((-3*c)^2 - 4*a^2))/2;

    % real and imaginary ranges
    real_part = linspace(low_real, up_real, acc);
    imag_part = linspace(low_imag, up_imag, acc);

    % meshgrid() call
    [Re, Im] = meshgrid(real_part, imag_part);

    % complex domain definition
    domain = Re + 1i*Im;

end
