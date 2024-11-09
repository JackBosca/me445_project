%{
main.m - Main script for the Aerodynamics ME-445 project. Selected paper: 
https://doi.org/10.1007/s00348-017-2429-4.

CALLED FUNCTIONS: src/

CALLED DATA FILES: -

REVISIONS:
- #v0, 04-11-2024, Release, Boscariol Jacopo

Changes: -
%}

clear
close all
clc

%% config call
run("config.m")

%% Optimal geometrical Joukowski transformation

% initial guess for parameters eta_origin, xi_origin, a
init_guess = [-0.1, 0.1, 1];

% error definition (to be specified in objective_geom)
% err_type_geom = 'sum_squared';

% minimize the objective function
opt_params_geom = fminsearch(@(params) objective_geom(params, x, y), ...
    init_guess);

% Joukowski function call with optimal params
[~, xj, yj] = joukowski_transform(opt_params_geom);

yj_interp_geom = profile_interpolator(xj, yj, x);

%% Optimal cl Joukowski transformation

% stall angles
alphas_p = 7;                           % (deg)
alphas_n = -8;                          % (deg)

% search linear interval in clW(alpha) curve
idx1 = find(fig8_clW(:, 1) > alphas_n, 1);
idx2 = find(fig8_clW(:, 1) > alphas_p, 1);
alphaW = deg2rad(fig8_clW(idx1:idx2, 1));
clW = fig8_clW(idx1:idx2, 2);

% search linear interval in clT(alpha) curve
idx1 = find(fig8_clT(:, 1) > alphas_n, 1);
idx2 = find(fig8_clT(:, 1) > alphas_p, 1);
alphaT = deg2rad(fig8_clT(idx1:idx2, 1));
clT = fig8_clT(idx1:idx2, 2);

% error definition (to be specified in objective_geom)
% err_type_cl = 'sum_squared';

% minimize the objective function
opt_params_clW = fminsearch(@(params) objective_cl(params, alphaW, clW), ...
    init_guess);
opt_params_clT = fminsearch(@(params) objective_cl(params, alphaT, clT), ...
    init_guess);

% Joukowski function call with optimal params
[~, xjW, yjW] = joukowski_transform(opt_params_clW);
[~, xjT, yjT] = joukowski_transform(opt_params_clT);

yj_interp_clW = profile_interpolator(xjW, yjW, x);
yj_interp_clT = profile_interpolator(xjT, yjT, x);

% cl potential flow function
cl = @(alpha, lambda, a) 2*pi*(alpha + lambda/a);

%% plots call
run("plots.m")
