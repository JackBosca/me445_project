function err = obj_cl_full(params, alpha, cl, err_type)
% obj_cl_full - Objective function for the Joukowski profile cl error 
% minimization problem, considering all params (eta, xi, a).
% 
% INPUTS: 
% - params(1), float: eta coordinate of the circle origin
% - params(2), float: xi coordinate of the circle origin
% - params(3), float: a parameter of the complex circle
% - alpha, float: alpha values to compute cl (rad)
% - cl, float: objective airfoil cl values
% - err_type, string: type of error to be minimized (default 'mean-abs')
%
% OUTPUTS:
% - err, float: error between Joukowski cl curve and objective one
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 10/11/24, Boscariol Jacopo
%               Changes: release.

    % Joukowski cl computation
    cl_j = 2*pi*(params(2)/(params(3) - params(1)) + alpha);

    if nargin < 4
        err_type = 'mean-abs';
    end

    % error computation
    switch err_type
        case 'sum-squared'
            err = sum((cl - cl_j).^2);
        case 'mean-squared'
            err = mean((cl - cl_j).^2);
        case 'mean-abs'
            err = mean(abs(cl - cl_j));
        otherwise
            error(['Error type should be either sum-squared, ' ...
                'mean-squared or mean-abs'])
    end

end
