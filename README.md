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

The `\main` folder contains the scripts that needs to be executed by the user and represents the logic of each theories (i.e. Potential flow and thin airfoil theory). Those scripts reproduces Figure 8 and 10 and are respectively :
- `main_fig8.m`
- `main_fig10.m`

These files calls the following functions that produces the graphs using basic Matlab plot function :
- `plots_fig8.m` 
- `plots_fig10.m` 

Both `main_fig8.m` and `main_fig10.m` calls functions sets in `\src` folder : 

- `config.m` : Loads external datas of the experimental results and NACA profil $$x,y$$ coordinates 
- `joukowski_transform.m` : Function performing the Joukowski transformation of a complex circle given its center $$(\chi_c,\eta_c)$$ and radius $$a$$ in the complex $$\zeta$$-plane
- `obj_geom.m` : Objective function for the Joukowski profile geometric error minimization problem (returns the error between the joukowski profile and the studied NACA profile)
- `profile_interpolator.m` : Function to interpolate the $$y_1$$ values of a certain airfoil into a grid $$x_2$$
- `obj_geom_eta.m` : Objective function for the Joukowski profile geometric error minimization problem with xi_origin, a being fixed (returns error between Joukowski profile and objective one)
- `get_L.m` : Function to ask the user values for the tail LE $$x$$-coordinate if he does not want to use the default ones (taken from the paper)
- `get_H.m` : Function to ask the user values for the tail LE $$y$$-coordinate if he does not want to use the default ones (taken from the paper)
- `get_domain` : function to automatically generate a complex domain suitable to perform the needed calculations and plots for the two-airfoils configuration. Upper and lower boundaries for the real and imaginary coordinates of the complex plane are computed by inverse trsnsforming the points at
   $$x_{\textrm{max}} = L + 3c,
   x_{\textrm{min}} = -3c,
   y_{\textrm{max}} = H + 3c,
   y_{\textrm{min}} = -3c$$
  in the Joukowski profile plane.
- `complex_vel.m` : function to compute the velocity field around the Joukowski complex circle and the corresponding Joukowski profile
- `vel_interpolator.m` : interpolates the velocity field $$(U, V)$$ in the point $$(L, H)$$ of the considered geometrical domani $$z_\textrm{domain}$$ (which does not need to be generated through meshgrid)
- `compute_coeffs.m` : function to compute $$A_0, A_1$$ coefficients of the thin aifoil expansion
- `obj_cl_full.m` : Objective function for the Joukowski profile $$C_l$$ error minimization problem, considering all params ($$\eta, \xi, a$$)
- `obj_cl_partial.m` : Objective function for the Joukowski profile $$C_l$$ error minimization problem, considering $$\xi$$ and $$a$$
- `cl_interpolator.m` : Function to interpolate the $$C_l$$ values of a certain airfoil into a new grid $$\alpha_2$$
