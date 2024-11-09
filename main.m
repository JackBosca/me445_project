%{
main.m - Main script for the Aerodynamics ME-445 project. Selected paper: 
https://doi.org/10.1007/s00348-017-2429-4.

CALLED FUNCTIONS: -

CALLED DATA FILES: -

REVISIONS:
- #v0, 04-11-2024, Release, Boscariol Jacopo

Changes: -
%}

clear
close all
clc

% config call
run("config.m")

% plot NACA 23012
figure
plot(x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'b')
axis equal
xlabel('x/c')
ylabel('y/c')
grid on

%% Searching an optimal geometrical Joukowski transformation

% initial guess for parameters eta_origin, xi_origin, a
init_guess = [-0.1, 0.1, 1];

% minimize the objective function
opt_params_geom = fminsearch(@(params) objective_geom(params, x, y), ...
    init_guess);

% Joukowski function call with optimal params
[~, xj, yj] = joukowski_transform(opt_params_geom);

yj_interp = profile_interpolator(xj, yj, x);

% add Joukowski plot to the NACA23012 one
hold on
plot(x, yj_interp, 'Marker', 'p')

%% Searching an optimal cl Joukowski transformation

% search linear interval in cl(alpha) curve
idx1 = find(fig8_clW(:, 1) > 0, 1);
idx2 = find(fig8_clW(:, 1) > 7, 1);
alpha = deg2rad(fig8_clW(idx1:idx2, 1));
cl = fig8_clW(idx1:idx2, 2);

% minimize the objective function
opt_params_cl = fminsearch(@(params) objective_cl(params, alpha, cl), ...
    init_guess);

% Joukowski function call with optimal params
[~, xj, yj] = joukowski_transform(opt_params_cl);

yj_interp = profile_interpolator(xj, yj, x);

% add Joukowski plot to the NACA23012 one
plot(x, yj_interp, 'Marker', 'p')

% dynamic legend
if ~exist('err_type_geom', 'var')
    err_type_geom = 'mean-abs';
end
lgd_text_geom = sprintf('Joukowski Airfoil - Geometry Opt by %s', ...
    err_type_geom);
lgd_text_cl = sprintf('Joukowski Airfoil - Cl Opt by %s', ...
    err_type_geom);
legend('NACA 23012', lgd_text_geom, lgd_text_cl)

%% wing cl

figure
plot(fig8_clW(:, 1), fig8_clW(:, 2), 'Marker', 'p')
grid on
xlabel('\alpha (deg)')

cl = @(alpha, lambda, a) 2*pi*(alpha + lambda/a);

alpha = deg2rad(0:1:10);

hold on
plot(rad2deg(alpha), cl(alpha, opt_params_geom(2), opt_params_geom(3)), ...
    'Marker', 'x')
plot(rad2deg(alpha), cl(alpha, opt_params_cl(2), opt_params_cl(3)), ...
    'Marker', 'x')

% dynamic legend
if ~exist('err_type_cl', 'var')
    err_type_cl = 'mean-abs';
end
lgd_text_geom = sprintf('Joukowski Airfoil - Geometry Opt by %s', ...
    err_type_cl);
lgd_text_cl = sprintf('Joukowski Airfoil - Cl Opt by %s', ...
    err_type_cl);
legend('NACA 23012 cl_W', lgd_text_geom, lgd_text_cl)
