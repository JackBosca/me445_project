function cl_interp = cl_interpolator(alpha1, cl1, alpha2)
% cl_interpolator - Function to interpolate the cl values of a certain
% airfoil into a new grid alpha2. 
% 
% INPUTS: 
% - alpha1, float: initial alphas grid (rad)
% - cl1, float: initial cl values of the airfoil
% - alpha2, float: final alphas grid (rad)
%
% OUTPUTS:
% - cl_interp, float: cl values of the airfoil on alpha2
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 05/11/24, Boscariol Jacopo
%               Changes: release.
    
    % interpolation
    cl_interp = [interp1(alpha1, cl1, alpha2, 'linear', 'extrap')];

end