function [A0, A1] = compute_coeffs(yc, c)
% compute_coeffs.m - function to compute A0, A1 coefficients of the thin 
% aifoil expansion.
% 
% INPUTS: 
% - yc, float: camber line coordinates
% - c, float: chord
%
% OUTPUTS:
% - A0, float: first thin airfoil expansion coefficient
% - A1, float: second thin airfoil expansion coefficient
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 26/11/24, Boscariol Jacopo
%               Changes: release.
    
    % defining t range
    t = linspace(0, pi, length(yc));
    
    % computing x based on the relation x/c = (1 - cos(t))/2
    x = c * (1 - cos(t)) / 2;
    
    % computing dyc_dx using central differences
    dyc_dx = gradient(yc) ./ gradient(x);
    
    % computing dx/dt = c*sin(t)/2
    dx_dt = c*sin(t)/2;
    
    % A0 integrand
    integrand = dyc_dx .* dx_dt;
    
    % trapz integration
    A0 = 1/pi*trapz(t, integrand);

    % A1 integrand
    integrand = dyc_dx .* dx_dt .* cos(t);
    
    % trapz integration
    A1 = 2/pi*trapz(t, integrand);

end
