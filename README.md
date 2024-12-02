# [ME-445] Aerodynamics - Group Project

## General Information 

This project is conducted as part of the EPFL Aerodynamics Master's course. The objective is to reproduce two figures from a selected external paper using two out of the three aerodynamic theories discussed in the lectures. 

In this project, Figure 8 and Figure 10 from the selected paper (https://doi.org/10.1007/s00348-017-2429-4), which examines the aerodynamic characteristics of a two-airfoil arrangement (a NACA 23012 wing and tail), will be reproduced using two different aerodynamic theories. Specifically:

- For Figure 8, the lift coefficient of a single wing will be computed using the <ins>**potential flow theory**</ins> :
<img src="/../main/data/fig8.jpg" width="300">

- For Figure 10, the effect of the wing on the tail's flow will be analyzed by computing the difference between the lift coefficients both with and without the presence of the wing using the <ins>**thin airfoil theory**</ins> :
<img src="/../main/data/fig10.jpg" width="300">

The data from these figures were extracted on October 29th, 2024, using the online tool WebPlotDigitizer (https://automeris.io/wpd/?v=5_2).

The authors of this project are:

- Aristide Paul Marc Vignon
- Haroun Naina
- Jacopo Boscariol

The entire project has been implemented in the MATLAB environment using its default functions, without reliance on any external modules.

## Data and file overview

Figure 8 and 10 are reproduced by running main_fig8.m and main_fig10.m, respectively. Each of these files represents the main logic of each theories (i.e. Potential flow and thin airfoil theory) and only calls the required functions to determine the wished quantities. All of the functions are set in the \src folder and are described below : 

- 
