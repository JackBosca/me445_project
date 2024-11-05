%{
main.m - Main script for the Aerodynamics ME-445 project. Selected paper: 
https://doi.org/10.1007/s00348-017-2429-4.

CALLED FUNCTIONS: -

CALLED DATA FILES: -

REVISIONS :
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

%% Searching a nice Joukowski Transformation

% initial guess for parameters eta_origin, xi_origin, a
init_guess = [-0.1, 0.1, 1];

% minimize the objective function
opt_params = fminsearch(@(params) objective_fct(params, x, y), ...
    init_guess);

% Joukowski function call with optimal params
[~, xj, yj] = joukowski_transform(opt_params);

yj_interp = profile_interpolator(xj, yj, x);

% add Joukowski plot to the NACA23012 one
hold on
plot(x, yj_interp, 'Marker', 'p')
legend('NACA 23012', 'Joukowski Airfoil')
