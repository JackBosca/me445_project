function [zeta_circle, profile] = joukowski_transform(eta_origin, xi_origin, a, acc)
% joukowski_transform - Function to perform the Joukowski transform of a
% complex circle given its origin and the a parameter.
% 
% INPUTS: 
% - eta_origin, float: eta coordinate of the circle origin
% - xi_origin, float: xi coordinate of the circle origin
% - a, float: a parameter of the complex circle
% - acc, int: accuracy for the circle discretization
%
% OUTPUTS:
% - zeta_circle, complex float: coordinates of the complex circle
% - profile, complex float: coordinates of the Joukowski profile
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 04/11/24, Boscariol Jacopo
%               Changes: release.

    % checking accuracy specification
    if nargin < 4
        acc = 101;
    end

    % defining circle center
    zeta_origin = eta_origin + 1i*xi_origin;
    
    % radii computation
    R = abs(zeta_origin - a);               % by vectorial definition of R
    R0 = abs(zeta_origin);
    delta = pi - angle(zeta_origin);
    
    % creating the complex circle
    zeta_circle = R*exp(1i*linspace(0, 2*pi, acc)) + R0*(-cos(delta) ...
        + 1i*sin(delta));
    
    % Joukowski transform
    profile = zeta_circle + a^2 ./ zeta_circle; 

end
