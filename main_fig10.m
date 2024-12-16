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
%close all
clc

%% config call
run("config.m")

%% User inputs and scalings

L = get_L(cW, L);
H = get_H(H);

% scaling NACA23012 x, y coordinates for wing
xW = cW*x;
yW = cW*y;

% scaling NACA23012 x, y coordinates for tail
xT = cT*x + L;
yT = cT*y + H;

% pre-processing wing coordinates for correct ycT computation
yT_ss = flip(yT(1:31));
xT_ss = flip(xT(1:31));

% interpolating pressure side coordinates on suction side ones
yT_ps = yT(31:end);
xT_ps = xT(31:end);

yT_ps = interp1(xT_ps, yT_ps, xT_ss, 'linear');
yT_ps(end) = yT(end);

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

%% Velocity field computation

% getting complex domain
domain = get_domain(opt_params(3), L, H, cT);

% getting cricle and profile velocity fields
[z_domain, U_c, V_c, U_p, V_p] = complex_vel(U_inf, alphaT, ...
    opt_params, domain);

% translating domain to have wing LE at (0, 0)
z_domain = z_domain + 2*opt_params(3);

% computing velocity magnitudes
vel_mag_c = sqrt(U_c.^2 + V_c.^2);
vel_mag_p = sqrt(U_p.^2 + V_p.^2);

% computing velocity directions
angles_c = rad2deg(atan2(V_c, U_c));
angles_p = rad2deg(atan2(V_p, U_p));

% retrieving velocity angles at LE of the tail
[U_LH, V_LH, alpha_effT] = vel_interpolator(z_domain, U_p, V_p, L, H);

%% Thin airfoil theory Cl computations

% computing camber line of the wing
ycT = (yT_ss + yT_ps)/2;

% compute_coeffs() call
[A0, A1] = compute_coeffs(ycT, cT);

% computing tail Cl at both single and 2-airfoils configurations
clT_S = 2*pi*(alphaT - A0) + pi*A1;
clT = 2*pi*(alpha_effT - A0) + pi*A1;

% computing fig 10 ratio
ratio = (clT - clT_S) ./ clT_S;

%% computing uncertainties on cl
[uxW, uxT, uyW, uyT] = RSS();

%% plots_fig10 call

% select alpha to use for velocity fields visualizations
idx = get_idx(alphaT);

run("plots_fig10.m")
