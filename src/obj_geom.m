function err = obj_geom(params, x, y, err_type)
% objective_geom - Objective function for the Joukowski profile geometric 
% error minimization problem.
% 
% INPUTS: 
% - params(1), float: eta coordinate of the circle origin
% - params(2), float: xi coordinate of the circle origin
% - params(3), float: a parameter of the complex circle
% - x, float: final grid for profile interpolator
% - y, float: objective airfoil coordinates
% - err_type, string: type of error to be minimized (default 'mean-abs')
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

    if nargin < 4
        err_type = 'mean-abs';
    end

    % error computation
    switch err_type
        case 'sum-squared'
            err = sum((y - yj_interp).^2);
        case 'mean-squared'
            err = mean((y - yj_interp).^2);
        case 'mean-abs'
            err = mean(abs(y - yj_interp));
        otherwise
            error(['Error type should be either sum-squared, ' ...
                'mean-squared or mean-abs'])
    end

end
