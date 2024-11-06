function err = objective_geom(params, x, y)
% objective_geom - Objective function for the Joukowski profile geometric 
% error minimization problem.
% 
% INPUTS: 
% - params(1), float: eta coordinate of the circle origin
% - params(2), float: xi coordinate of the circle origin
% - params(3), float: a parameter of the complex circle
% - x, float: final grid for profile interpolator
% - y, float: objective airfoil coordinates
%
% OUTPUTS:
% - err, float: error between Joukowski profile and objective one
%
% CALLED FUNCTIONS: joukowski_transform, profile_interpolator
%
% REVISIONS:
% - #v0 05/11/24, Boscariol Jacopo
%               Changes: release.

    % Joukowski trasnformation
    [~, xj, yj] = joukowski_transform(params);

    % interpolation on x grid
    yj_interp = profile_interpolator(xj, yj, x);

    % error computation
    % err = sum((y - yj_interp).^2);
    err = mean(abs(y - yj_interp));

end
