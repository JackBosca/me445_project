%{
plots.m - Plots script for the Aerodynamics ME-445 project. Selected paper: 
https://doi.org/10.1007/s00348-017-2429-4.

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

%% Airfoils
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
plot(x, yj_interp_clW, 'Marker', 'p', 'Color', colW)
plot(x, yj_interp_clT, 'Marker', 'p', 'Color', colT)

% dynamic legend
if ~exist('err_type_geom', 'var')
    err_type_geom = 'mean-abs';
end
lgd_text_geom = sprintf('Joukowski Airfoil - Geometry Opt by %s', ...
    err_type_geom);
lgd_text_clW = sprintf('Joukowski Wing Airfoil - cl Opt by %s', ...
    err_type_geom);
lgd_text_clT = sprintf('Joukowski Tail Airfoil - cl Opt by %s', ...
    err_type_geom);

legend('NACA 23012', lgd_text_geom, lgd_text_clW, lgd_text_clT)

%% cl curves

figure
plot(fig8_clW(:, 1), fig8_clW(:, 2), 'Marker', 'p', 'Color', 'black')
hold on
plot(fig8_clT(:, 1), fig8_clT(:, 2), 'Marker', 'square', 'Color', 'black')
grid on
xlabel('\alpha (deg)')

plot(rad2deg(alphaW), cl(alphaW, opt_params_geom(2), opt_params_geom(3)), ...
    'Marker', 'x', 'Color', col_geom)
plot(rad2deg(alphaW), cl(alphaW, opt_params_clW(2), opt_params_clW(3)), ...
    'Marker', 'x', 'Color', colW)
plot(rad2deg(alphaT), cl(alphaT, opt_params_clT(2), opt_params_clT(3)), ...
    'Marker', 'x', 'Color', colT)

% dynamic legend
if ~exist('err_type_cl', 'var')
    err_type_cl = 'mean-abs';
end
lgd_text_geom = sprintf('Joukowski Airfoil - Geometry Opt by %s', ...
    err_type_cl);
lgd_text_clW = sprintf('Joukowski Wing Airfoil - cl Opt by %s', ...
    err_type_geom);
lgd_text_clT = sprintf('Joukowski Tail Airfoil - cl Opt by %s', ...
    err_type_geom);

legend('NACA 23012 cl_W', 'NACA 23012 cl_T', lgd_text_geom, ...
    lgd_text_clW, lgd_text_clT)
title('Fig. 8')
