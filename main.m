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
init_guess_geom = [0.1, 0.1, 1];

% error definition
err_type_geom = 'mean-squared';

% minimize the objective function
opt_params_geom = fminsearch(@(params) obj_geom(params, x, y, ... 
    err_type_geom), init_guess_geom);

% Joukowski function call with optimal params
[~, xj, yj] = joukowski_transform(opt_params_geom);

yj_interp_geom = profile_interpolator(xj, yj, x);

%% Optimal partial (xi, a) cl Joukowski transformation

% stall angles
alphas_pW = 7;                           % (deg)
alphas_nW = -8;                          % (deg)
alphas_pT = 5;                           % (deg)
alphas_nT = -8;                          % (deg)

% search linear interval in clW(alpha) curve
idx1 = find(fig8_clW(:, 1) > alphas_nW, 1);
idx2 = find(fig8_clW(:, 1) > alphas_pW, 1);
alphaW = deg2rad(fig8_clW(idx1:idx2, 1));
clW = fig8_clW(idx1:idx2, 2);

% search linear interval in clT(alpha) curve
idx1 = find(fig8_clT(:, 1) > alphas_nT, 1);
idx2 = find(fig8_clT(:, 1) > alphas_pT, 1);
alphaT = deg2rad(fig8_clT(idx1:idx2, 1));
clT = fig8_clT(idx1:idx2, 2);

% initial guess
init_guess_cl = [0.1, 0.1, 1];

% error definition
err_type_cl = 'sum-squared';

% minimize the objective function
opt_params_clW_partial = fminsearch(@(params) obj_cl_partial(params, alphaW, clW, ...
    err_type_cl), init_guess_cl(2:3));
opt_params_clT_partial = fminsearch(@(params) obj_cl_partial(params, alphaT, clT, ...
    err_type_cl), init_guess_cl(2:3));

% optimal xi, a have been found, now search for best eta to find best 
% geometrical fit of the NACA profile
etaW = fminsearch(@(eta) obj_geom_eta(eta, opt_params_clW_partial(1), ...
    opt_params_clW_partial(2), x, y, err_type_geom), init_guess_cl(1));
etaT = fminsearch(@(eta) obj_geom_eta(eta, opt_params_clT_partial(1), ...
    opt_params_clT_partial(2), x, y, err_type_geom), init_guess_cl(1));

opt_params_clW_partial = [etaW, opt_params_clW_partial];
opt_params_clT_partial = [etaT, opt_params_clT_partial];

% Joukowski function call with optimal params
[~, xjW, yjW] = joukowski_transform(opt_params_clW_partial);
[~, xjT, yjT] = joukowski_transform(opt_params_clT_partial);

yj_interp_clW_partial = profile_interpolator(xjW, yjW, x);
yj_interp_clT_partial = profile_interpolator(xjT, yjT, x);

% cl potential flow function (apporx lambda/a --> 0 is valid for the
% NACA 23012)
cl_partial = @(alpha, lambda, a) 2*pi*(alpha + lambda/a);

%% Optimal full (eta, xi, a) cl Joukowski transformation

% initial guess
init_guess_cl = [0.1, 0.1, 1];

% error definition
err_type_cl = 'sum-squared';

% minimize the objective function
opt_params_clW_full = fminsearch(@(params) obj_cl_full(params, alphaW, clW, ...
    err_type_cl), init_guess_cl);
opt_params_clT_full = fminsearch(@(params) obj_cl_full(params, alphaT, clT, ...
    err_type_cl), init_guess_cl);

% Joukowski function call with optimal params
[~, xjW, yjW] = joukowski_transform(opt_params_clW_full);
[~, xjT, yjT] = joukowski_transform(opt_params_clT_full);

yj_interp_clW_full = profile_interpolator(xjW, yjW, x);
yj_interp_clT_full = profile_interpolator(xjT, yjT, x);

% cl potential flow function (apporx lambda/a --> 0 is valid for the
% NACA 23012)
cl_full = @(alpha, params) 2*pi*sin(alpha + ...
    params(2)/(params(3) - params(1)));

%% plots call
run("plots.m")
