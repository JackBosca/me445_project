function [zeta_circle, xj, yj] = joukowski_transform(params, acc)
% joukowski_transform - Function to perform the Joukowski transform of a
% complex circle given its origin and the a parameter.
% 
% INPUTS: 
% - params(1), float: eta coordinate of the circle origin
% - params(2), float: xi coordinate of the circle origin
% - params(3), float: a parameter of the complex circle
% - acc, int: accuracy for the circle discretization
%
% OUTPUTS:
% - zeta_circle, complex float: coordinates of the complex circle
% - xj, float: x coordinates of the Joukowski profile
% - yj, float: y coordinates of the Joukowski profile
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

    % extracting params
    eta_origin = params(1);
    xi_origin = params(2);
    a = params(3);

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

    % traslating profile to have LE at (0, 0)
    profile = profile - min(real(profile));

    % dividing by profile chord
    xj = real(profile)/(max(real(profile)) - min(real(profile)));
    yj = imag(profile)/(max(real(profile)) - min(real(profile)));

end
