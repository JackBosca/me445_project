%{
config.m - Configuration script for the Aerodynamics ME-445 project. 

CALLED FUNCTIONS: -

CALLED DATA FILES: data/

REVISIONS:
- #v0, 04-11-2024, Release, Boscariol Jacopo

Changes: -
%}

% paths
data_path = "data/";
addpath(genpath(data_path))
src_path = "src/";
addpath(genpath(src_path))

% load NACA 23012 data
load("x.mat")                           % x/c
load("y.mat")                           % y/c

% load figures data
load("fig8_clW.mat")
load("fig8_clT.mat")
load("fig10_T.mat")

% geometric and aerodynamic data
c = 0.05;                               % (m)
cW = c;                                 % wing chord
cT = c/2;                               % tail chord
L = 3*c;                                % longitudinal W-T distance
H = c;                                  % vertical W-T distance
U_inf = [17.5, 35];                     % (m/s), [2-airfoils, wing]
Re_c = [5.83e4, 1.16e5];                % (-), chord Re for wing [17.5, 35] 

% getting kinematic viscosity
v = U_inf*cW ./ Re_c;                   % (m^2/s)
