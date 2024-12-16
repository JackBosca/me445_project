******** [ME-445] Aerodynamics - Project (Group 24) ********

**** General Information ****

This project is conducted as part of the EPFL Aerodynamics Master's course. The objective is to reproduce two figures from a selected external paper using two out of the three aerodynamic theories discussed in the lectures.

In this project, Figure 8 and Figure 10 from the selected paper (https://doi.org/10.1007/s00348-017-2429-4), which examines the aerodynamic characteristics of a two-airfoil arrangement (a NACA 23012 wing and tail), will be reproduced using two different aerodynamic theories. Specifically:

- For Figure 8, the lift coefficient of a single wing will be computed using the potential flow theory.
- For Figure 10, the effect of the wing on the tail's flow will be analyzed by computing the difference between the lift coefficients both with and without the presence of the wing using the thin airfoil theory.

The authors of this project are:
- Aristide Paul Marc Vignon;
- Haroun Naina;
- Jacopo Boscariol.

The entire project has been implemented in the MATLAB environment using its default functions, without reliance on any external modules.





**** Data and File Overview ****

The project is run by the execution of two main scripts. These scripts reproduce Figure 8 and Figure 10 by applying the corresponding theories (Potential flow and thin airfoil theory, respectively) and are located outside any folder:
- main_fig8.m;
- main_fig10.m.

These files call the configuration script config.m and the following functions that produce the graphs using basic MATLAB plot functions:
- plots_fig8.m;
- plots_fig10.m.

Both main_fig8.m and main_fig10.m call functions set in the \src folder:
- joukowski_transform.m: Function performing the Joukowski transformation of a complex circle given its center (chi_c, eta_c) and radius a in the complex zeta-plane;
- obj_geom.m: Objective function for the Joukowski profile geometric error minimization problem (returns the error between the Joukowski profile and the studied NACA profile);
- profile_interpolator.m: Function to interpolate the y_1 values of a certain airfoil into a grid x_2;
- obj_geom_eta.m: Objective function for the Joukowski profile geometric error minimization problem with xi-origin and a fixed (returns error between the Joukowski profile and the objective one);
- get_L.m: Function to ask the user values for the tail LE x-coordinate if they do not want to use the default ones (taken from the paper);
- get_H.m: Function to ask the user values for the tail LE y-coordinate if they do not want to use the default ones (taken from the paper);
- get_domain: Function to automatically generate a complex domain suitable to perform the needed calculations and plots for the two-airfoil configuration. Upper and lower boundaries for the real and imaginary coordinates of the complex plane are computed by inverse transforming points in the Joukowski profile plane;
- complex_vel.m: Function to compute the velocity field around the Joukowski complex circle and the corresponding Joukowski profile;
- vel_interpolator.m: Interpolates the velocity field (U, V) in the point (L, H) of the considered geometrical domain z_domain (which does not need to be generated through meshgrid);
- compute_coeffs.m: Function to compute A_0, A_1 coefficients of the thin airfoil expansion;
- obj_cl_full.m: Objective function for the Joukowski profile Cl error minimization problem, considering all parameters (eta, xi, a);
- obj_cl_partial.m: Objective function for the Joukowski profile Cl error minimization problem, considering xi and a;
- cl_interpolator.m: Function to interpolate the Cl values of a certain airfoil into a new grid alpha_2.

The \data folder contains the paper's experimental results and airfoil x, y coordinates:
- NACA23012.txt: x, y coordinates of the studied airfoil (2x61);
  - x.mat: x coordinate (1x61);
  - y.mat: y coordinate (1x61);
- fig8_clT.mat: Cl(alpha) of tail extracted from Figure 8 (2x80);
- fig8_clW.mat: Cl(alpha) of wing extracted from Figure 8 (2x80);
- fig9_clT.mat: Cl(alpha) of tail extracted from Figure 9 (2x98);
- fig10_T.mat: (Cl - Cl,s)/Cl,s as a function of alpha (2x51);
- fig8.jpg: Figure 8;
- fig10.jpg: Figure 10;
- arrangement.jpg: Figure representing the arrangement of both airfoils.






**** Methodological Information ****

The user only needs to run the following scripts, everything else is completely automatic:
- main_fig8.m;
- main_fig10.m.

The execution of script main_fig8.m does not ask any input from the user and computes coefficient Cl according to the paper's experimental setup, using potential flow theory.

- Input:
  - None;
- Output:
  - Airfoil partial optimization;
  - Airfoil full optimization;
  - Figure comparing Cl(alpha): Experimental, partial optimization;
  - Figure comparing Cl(alpha): Experimental, full optimization.

The execution of script main_fig10.m asks the user the values of L and H. These variables define the position of the tail as described in the paper. It also asks a value of undisturbed α to be used for plots. The script then computes (Cl - Cl,s)/Cl,s for the specified position using thin airfoil theory.

- Input:
  - L Value;
  - H Value;
  - α Value;
- Output:
  - Velocity streamlines and magnitude contour around the circle in zeta-plane;
  - Velocity streamlines and direction contour around the circle in zeta-plane;
  - Velocity magnitude contour around Joukowski profile;
  - Velocity direction contour around Joukowski profile;
  - Figure comparing (Cl - Cl,s)/Cl,s as a function of alpha: Experimental (Paper fig10), TA Theory;
  - Figure comparing Cl as a function of alpha: Experimental (Paper fig9), TA Theory.
