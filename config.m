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
load("fig9_clT.mat")
load("fig10_T.mat")

% geometric and aerodynamic data
c = 0.05;                               % (m)
cW = c;                                 % wing chord
cT = c/2;                               % tail chord
L = 3*c;                                % longitudinal W-T distance
H = c;                                  % vertical W-T distance
U_inf = 17.5;                     % (m/s)
Re_c = 5.83e4;                % (-) chord Re for wing 

% getting kinematic viscosity
v = U_inf*cW / Re_c;                   % (m^2/s)

% stall angles
alphas_pW = 7;                           % (deg)
alphas_nW = -8;                          % (deg)
alphas_pT = 5;                           % (deg)
alphas_nT = -8;                          % (deg)
