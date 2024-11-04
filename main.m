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
title('NACA 23012 Airfoil Profile')
grid on

%% Searching a Joukowski Transformation

% defining circle center
eta_origin = -0.1; 
xi_origin = 0.1; 

% defining a parameter
a = 1.2;

% Joukowski function call
[zeta_circle, xj, yj] = joukowski_transform(eta_origin, xi_origin, a);

% interpolating yj onto x
idx_j = find(xj == 0, 1);
idx_x = find(x == 0, 1);

yj_interp1 = interp1(xj(1:idx_j), yj(1:idx_j), x(1:idx_x), ...
    'linear', 'extrap');
yj_interp2 = interp1(xj((idx_j + 1):end), yj((idx_j + 1):end), ...
    x((idx_x + 1):end), 'linear', 'extrap');

yj_interp = [yj_interp1, yj_interp2];

% add Joukowski plot to the NACA23012 one
hold on
plot(x, yj_interp, 'Marker', 'p')
