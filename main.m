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
[zeta_circle, profile] = joukowski_transform(eta_origin, xi_origin, a);

% plots
figure 

subplot(1, 2, 1)
plot(complex(0), 'k*')                          % (0, 0)
hold on
plot(complex(eta_origin + 1i*xi_origin), 'bo')  % center of the circle 
plot(zeta_circle, 'm')
axis equal tight
axis([-1.5, 1.5, -1.5, 1.5]) 
grid on
xlabel('\eta')
ylabel('\xi')

subplot(1,2,2)
plot(profile, 'm')
hold on
axis equal tight
axis([-3, 3, -3, 3]) 
grid on
xlabel('x')
ylabel('y') 
