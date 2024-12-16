%{
plots_fig8.m - Plots script for the Aerodynamics ME-445 project, figure 8. 
Selected paper: https://doi.org/10.1007/s00348-017-2429-4.

CALLED FUNCTIONS: -

CALLED DATA FILES: -

REVISIONS:
- #v0, 09-11-2024, Release, Boscariol Jacopo

Changes: -
%}

% set colors
col_geom = 'blue';
colW = 'red';
colT = 'green';

%% Airfoils (partial cl optimization)
% plot NACA 23012
figure
plot(x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 3, ...
    'MarkerFaceColor', 'black', 'Color', 'black')
axis equal
xlabel('x/c')
ylabel('y/c')
grid on

% add opt geometrical Joukowski plot
hold on
plot(x, yj_interp_geom, 'Marker', 'p', 'Color', col_geom)

% add opt cl Joukowski plot
plot(x, yj_interp_clW_partial, 'Marker', 'p', 'Color', colW)
plot(x, yj_interp_clT_partial, 'Marker', 'p', 'Color', colT)

% dynamic legend
if ~exist('err_type_geom', 'var')
    err_type_geom = 'mean-abs';
end
lgd_text_geom = sprintf('Joukowski Airfoil - Geometry Opt by %s', ...
    err_type_geom);
lgd_text_clW = sprintf('Joukowski Wing Airfoil - c_l Opt by %s', ...
    err_type_cl);
lgd_text_clT = sprintf('Joukowski Tail Airfoil - c_l Opt by %s', ...
    err_type_cl);

legend('NACA 23012', lgd_text_geom, lgd_text_clW, lgd_text_clT)
title('Partial c_l Optimization')

%% cl (partial optimization) curves

figure
plot(fig8_clW(:, 1), fig8_clW(:, 2), 'k^')
hold on
plot(fig8_clT(:, 1), fig8_clT(:, 2), 'ksquare')
grid on
xlabel('\alpha (deg)')
ylabel('c_l (-)')

plot(rad2deg(alphaW), cl_partial(alphaW, opt_params_geom(2), ...
    opt_params_geom(3)), 'Marker', 'p', 'Color', col_geom)
plot(rad2deg(alphaW), cl_partial(alphaW, opt_params_clW_partial(2), ...
    opt_params_clW_partial(3)), 'Marker', 'p', 'Color', colW)
plot(rad2deg(alphaT), cl_partial(alphaT, opt_params_clT_partial(2), ...
    opt_params_clT_partial(3)), ...
    'Marker', 'p', 'Color', colT)

% dynamic legend
if ~exist('err_type_cl', 'var')
    err_type_cl = 'mean-abs';
end
lgd_text_geom = sprintf('Joukowski Airfoil - Geometry Opt by %s', ...
    err_type_geom);
lgd_text_clW = sprintf('Joukowski Wing Airfoil - c_l Opt by %s', ...
    err_type_cl);
lgd_text_clT = sprintf('Joukowski Tail Airfoil - c_l Opt by %s', ...
    err_type_cl);

legend('NACA 23012 c_l_W', 'NACA 23012 c_l_T', lgd_text_geom, ...
    lgd_text_clW, lgd_text_clT)
title('Fig. 8, Partial c_l Optimization')

%% Airfoils (full cl optimization)
% plot NACA 23012
figure
plot(x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 3, ...
    'MarkerFaceColor', 'black', 'Color', 'black')
axis equal
xlabel('x/c')
ylabel('y/c')
grid on

% add opt geometrical Joukowski plot
hold on
plot(x, yj_interp_geom, 'Marker', 'p', 'Color', col_geom)

% add opt cl Joukowski plot
plot(x, yj_interp_clW_full, 'Marker', 'p', 'Color', colW)
plot(x, yj_interp_clT_full, 'Marker', 'p', 'Color', colT)

% dynamic legend
if ~exist('err_type_geom', 'var')
    err_type_geom = 'mean-abs';
end
lgd_text_geom = sprintf('Joukowski Airfoil - Geometry Opt by %s', ...
    err_type_geom);
lgd_text_clW = sprintf('Joukowski Wing Airfoil - c_l Opt by %s', ...
    err_type_cl);
lgd_text_clT = sprintf('Joukowski Tail Airfoil - c_l Opt by %s', ...
    err_type_cl);

legend('NACA 23012', lgd_text_geom, lgd_text_clW, lgd_text_clT)
title('Full c_l Optimization')

%% cl (full optimization) curves

figure
plot(fig8_clW(:, 1), fig8_clW(:, 2), 'k^')
hold on
plot(fig8_clT(:, 1), fig8_clT(:, 2), 'ksquare')
grid on
xlabel('\alpha (deg)')
ylabel('c_l (-)')

% geometry opt cl(alpha)
plot(rad2deg(alphaW), cl_full(alphaW, opt_params_geom), ...
    'Marker', 'o', 'MarkerSize', 3, 'Color', col_geom)

% wing cl opt cl(alpha)
plot(rad2deg(alphaW), cl_full(alphaW, opt_params_clW_full), ...
    'Marker', 'o', 'MarkerSize', 3, 'Color', colW)

% tail cl opt cl(alpha)
plot(rad2deg(alphaT), cl_full(alphaT, opt_params_clT_full), ...
    'Marker', 'o', 'MarkerSize', 3, 'Color', colT)

% adding errorbars
errorbar(rad2deg(alphaW), cl_full(alphaW, opt_params_geom), ...
    uyW, uyW, rad2deg(uxW), rad2deg(uxW), 'Color', col_geom, 'CapSize', 3)
errorbar(rad2deg(alphaW), cl_full(alphaW, opt_params_clW_full), ...
    uyW, uyW, rad2deg(uxW), rad2deg(uxW), 'Color', colW, 'CapSize', 3)
errorbar(rad2deg(alphaT), cl_full(alphaT, opt_params_clT_full), ...
    uyT, uyT, rad2deg(uxT), rad2deg(uxT), 'Color', colT, 'CapSize', 3)

% dynamic legend
if ~exist('err_type_cl', 'var')
    err_type_cl = 'mean-abs';
end
lgd_text_geom = sprintf('Joukowski Airfoil - Geometry Opt by %s', ...
    err_type_geom);
lgd_text_clW = sprintf('Joukowski Wing Airfoil - c_l Opt by %s', ...
    err_type_cl);
lgd_text_clT = sprintf('Joukowski Tail Airfoil - c_l Opt by %s', ...
    err_type_cl);

legend('NACA 23012 c_l_W', 'NACA 23012 c_l_T', lgd_text_geom, ...
    lgd_text_clW, lgd_text_clT)
title('Fig. 8, Full c_l Optimization')
