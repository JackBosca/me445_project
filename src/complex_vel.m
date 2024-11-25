function [z_domain, u_c, v_c, u_p, v_p] = complex_vel(U_inf, alpha, params, domain)
% complex_vel - function to compute the velocity field around the Joukowski
% complex circle and the corresponding Joukowski profile.
% 
% INPUTS: 
% - U_inf, float: asymptotic velocity
% - alpha, float: asymptotic angle of attack
% - params(1), float: eta coordinate of the circle origin
% - params(2), float: xi coordinate of the circle origin
% - params(3), float: a parameter of the complex circle
% - domain, complex float: domain in the complex circle plane
%
% OUTPUTS:
% - z_domain, complex float: domain in the Joukowski profile plane
% - u_c, float: x velocity in the circle plane
% - v_c, float: y velocity in the circle plane
% - u_p, float: x velocity in the profile plane
% - v_p, float: y velocity in the profile plane
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 25/11/24, Boscariol Jacopo
%               Changes: release.
    
    % params unpacking
    eta = params(1);
    xi = params(2);
    a = params(3);

    % retrieving R
    R = sqrt((a - eta)^2 + xi^2);

    % mask to exclude points inside the complex circle
    R2 = (a - eta)^2 + xi^2;
    dist2 = (real(domain) - eta).^2 + (imag(domain) - xi).^2;
    mask = dist2 < R2;

    % retrieving gamma from Kutta condition
    gamma = -1i*2*pi*(a - eta -1i*xi)*(U_inf*exp(-1i*alpha) - ...
        U_inf*R^2*exp(1i*alpha)/(a - eta - 1i*xi)^2);

    % complex potential derivative
    dW_dz = @(zeta) U_inf*exp(-1i*alpha) - ...
        U_inf*R^2*exp(1i*alpha) ./ (zeta - eta - 1i*xi).^2 - ...
        1i*gamma ./ (2*pi*(zeta - eta - 1i*xi));

    % inverse derivative
    dz_dzeta = @(zeta) 1 - a^2 ./ (zeta.^2);

    % transforming in physical domain
    z_domain = domain + a^2 ./ domain;

    % retrieving circle velocities
    u_c = real(dW_dz(domain));
    v_c = -imag(dW_dz(domain));

    % filtering points of the domain
    u_c(mask) = NaN;
    v_c(mask) = NaN;

    % retrieving profile velocities
    u_p = real(dW_dz(z_domain) ./ dz_dzeta(z_domain));
    v_p = -imag(dW_dz(z_domain) ./ dz_dzeta(z_domain));

    % filtering points of the domain
    u_p(mask) = NaN;
    v_p(mask) = NaN;

end
