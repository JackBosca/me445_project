function y_interp = profile_interpolator(x1, y1, x2)
% profile_interpolator - Function to interpolate the y1 values of a certain
% airfoil into a grid x2.
% 
% INPUTS: 
% - x1, float: initial grid
% - y1, float: initial y coordinates of the airfoil
% - x2, float: final grid
%
% OUTPUTS:
% - y_interp, float: y coordinates of the airfoil on x2
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 05/11/24, Boscariol Jacopo
%               Changes: release.

    % indexes to split x grids
    idx_1 = find(x1 == 0, 1);
    idx_2 = find(x2 == 0, 1);
    
    % interpolation
    y_interp = [interp1(x1(1:idx_1), y1(1:idx_1), x2(1:idx_2), ...
        'linear', 'extrap'), ...
        interp1(x1((idx_1 + 1):end), y1((idx_1 + 1):end), ...
        x2((idx_2 + 1):end), 'linear', 'extrap')];

end