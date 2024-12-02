# [ME-445] Aerodynamics - Project (Group 24)

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

The project is run by the execution of two main script. Those scripts represents the logic of each theories (i.e. Potential flow and thin airfoil theory) and are repertioried in the `\main` folder :
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

The `\data` folder is, as its name indicates, the papers experimental result and airfoil $$x,y$$ coordinates : 
- `NACA23012.txt` : $$x,y$$ coordinates of the studied airfoil (2x61)
   - `x.mat` : $$x$$ coordinate (1x61)
   - `y.mat` : $$y$$ coordinate (1x61)
- `fig8_clT.mat` : $$C_l(\alpha)$$ of tail extracted from figure 8 (2x80)
- `fig8_clW.mat` : $$C_l(\alpha)$$ of wing extracted from figure 8 (2x80)
- `fig10_T.mat` : $$(C_l - C_{l,s})/C_{l,s}$$ as a funnction of $$\alpha$$ (2x51)
   - `fig10_tail_CL.txt` : $$(C_l - C_{l,s})/C_{l,s}$$ of tail from figure 10 (1x51)
   - `fig10_tail_alpha.txt` : $$\alpha$$ of tail from figure 10 (1x51)
- `fig8.jpg` : Figure 8
- `fig10.jpg` : Figure 10
- `arrangement.jpg` : Figure representing the arrangement of both airfoil

The `\_old` folder contains previous work that led to the establishing of the final scripts.

## Methodological information

- The execution of script `main_fig8.m` does not ask any input from the user and computes coefficient $$C_l$$ according to the paper experimental setup, using potential flow theory.
   - <ins>**Input**</ins> :
      - [-]
   - <ins>**Output**</ins> :
      - Airfoil partial optimization
      - Airfoil full optimization
      - Figure comparing $$C_l(\alpha)$$ : Experimental, partial optimization
      - Figure comparing $$C_l(\alpha)$$ : Experimental, full optimization
- The execution of script `main_fig10.m` asks the user the values of $$L$$ and $$H$$. Those variables defines the position of the tail as described in the paper (see figure below). The script, then, computes $$(C_l - C_{l,s})/C_{l,s}$$ for the specified position using thin airfoil theory.
   - <ins>**Input**</ins> :
      - $$L$$ Value
      - $$H$$ Value
   - <ins>**Output**</ins> :
      - Velocity streamlines and magnitude contour around circle in $$\zeta$$-plane
      - Velocity streamlines and direction contour around circle in $$\zeta$$-plane 
      - Velocity magnitude contour around Joukowski profile
      - Velocity direction contour around Joukowski profile
      - Figure comparing $$(C_l - C_{l,s})/C_{l,s}$$ as a function of $$\alpha$$ : Experimental, TA Theory (for the specified location of the airfoil)

<img src="/../main/data/arrangement.jpg" width="300">
