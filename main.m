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
