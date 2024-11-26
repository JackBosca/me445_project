%{
main_fig10.m - Main script for the Aerodynamics ME-445 project, figure 10. 
Selected paper: https://doi.org/10.1007/s00348-017-2429-4.

CALLED FUNCTIONS: src/

CALLED DATA FILES: -

REVISIONS:
- #v0, 26-11-2024, Release, Boscariol Jacopo

Changes: -
%}

clear
close all
clc

%% config call
run("config.m")

% scaling NACA23012 x, y coordinates for wing
xW = cW*x;
yW = cW*y;

% scaling NACA23012 x, y coordinates for tail
xT = cT*x + L;
yT = cT*y + H;

% getting considered angles of attack
alphaT = deg2rad(fig10_T(:, 1));

%% Optimal geometrical Joukowski transformation

% initial guess for parameters eta_origin, xi_origin, a
init_guess = [0.1, 0.1, 1];

% error definition
err_type = 'mean-squared';

% minimize the objective function
opt_params = fminsearch(@(params) obj_geom(params, x, y, ... 
    err_type), init_guess);

% scaling parameters to have correct wing chord
opt_params = opt_params * cW / (4*opt_params(3));

% Joukowski function call with optimal params
[zeta_circle, xj, yj] = joukowski_transform(opt_params);

%% 
mR = -0.03;
MR = 0.2;
mI = -0.03;
MI = 0.08;

realPart = linspace(mR, MR, 1000);
imagPart = linspace(mI, MI, 1000);

[Re, Im] = meshgrid(realPart, imagPart);
domain = Re + 1i*Im;

[z_domain, U_c, V_c, U_p, V_p] = complex_vel(U_inf(1), alphaT, ...
    opt_params, domain);

velocity_mag_c = sqrt(U_c.^2 + V_c.^2);
velocity_mag_p = sqrt(U_p.^2 + V_p.^2);

angles_c = atan2(V_c, U_c);
angles_p = atan2(V_p, U_p);

angles_c = rad2deg(angles_c);
angles_p = rad2deg(angles_p);

%% CIRCLE
figure
plot(real(zeta_circle), imag(zeta_circle), 'Color', 'k')
axis equal
hold on

% adding vel magnitude contours
contourf(real(domain), imag(domain), velocity_mag_c(:, :, 1), 20, 'LineColor', 'none');
colorbar;

% adding streamlines with black color and finer density
h_stream = streamslice(real(domain), imag(domain), U_c(:, :, 1), V_c(:, :, 1), 3);
set(h_stream, 'Color', 'k'); 

title('Velocity Magnitude Around Joukowski Circle');
xlabel('Re(z)');
ylabel('Im(z)');

figure
plot(real(zeta_circle), imag(zeta_circle), 'Color', 'k')
axis equal
hold on

% adding vel angle contours
contourf(real(domain), imag(domain), angles_c(:, :, 1), 20, 'LineColor', 'none');
colorbar;

% adding streamlines with black color and finer density
h_stream = streamslice(real(domain), imag(domain), U_c(:, :, 1), V_c(:, :, 1), 3);
set(h_stream, 'Color', 'k'); 

title('Velocity Angle Around Joukowski Circle');
xlabel('Re(z)');
ylabel('Im(z)');

%% PROFILE

% scaling z_domain to have a unitary chord profile
% z_domain = z_domain/(4*opt_params_geom(3));

% translating domain to have wing LE at (0, 0)
z_domain = z_domain + 2*opt_params(3);

figure;
contourf(real(z_domain), imag(z_domain), velocity_mag_p(:, :, 1), 20, 'LineColor', 'none');
colorbar;
hold on;
plot(xW, yW, '-o', 'LineWidth', 1.5, 'MarkerSize', 3, ...
    'MarkerFaceColor', 'black', 'Color', 'black')
title('Velocity Magnitude Around Joukowski Profile');
xlabel('Re(z)');
ylabel('Im(z)');
axis equal;

figure;
contourf(real(z_domain), imag(z_domain), angles_p(:, :, 1), 20, 'LineColor', 'none');
colorbar;
hold on;
plot(xW, yW, '-o', 'LineWidth', 1.5, 'MarkerSize', 3, ...
    'MarkerFaceColor', 'black', 'Color', 'black')
plot(xT, yT, '-o', 'LineWidth', 1.5, 'MarkerSize', 3, ...
    'MarkerFaceColor', 'black', 'Color', 'black')
title('Velocity Angle Around Joukowski Profile');
xlabel('Re(z)');
ylabel('Im(z)');
axis equal;

%% compute camber-line
ycW = yW/2;

[A0, A1] = compute_coeffs(ycW, cW);

clT_S = 2*pi*(alphaT - A0) + pi*A1;

figure
plot(fig8_clW(:, 1), fig8_clW(:, 2), 'Marker', '^', 'Color', 'black')
hold on
plot(fig8_clT(:, 1), fig8_clT(:, 2), 'Marker', 'square', 'Color', 'black')
grid on
xlabel('\alpha (deg)')
ylabel('c_l (-)')

% geometry opt cl(alpha)
plot(rad2deg(alphaT), clT_S, ...
    'Marker', 'o', 'MarkerSize', 3)

%% velocity at (L, H)

[U_LH, V_LH, angles] = vel_interpolator(z_domain, U_p, V_p, L, H);

%%
clT = 2*pi*(angles - A0) + pi*A1;

q = (clT - clT_S) ./ clT_S;

idx = true(1, length(alphaT));
idx(21) = false;

figure
plot(fig10_T(:, 1), fig10_T(:, 2), 'ko')
hold on
grid on
plot(rad2deg(alphaT(idx)), q(idx))
