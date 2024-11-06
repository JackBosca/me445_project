function err = objective_cl(params, alpha, cl)
% objective_cl - Objective function for the Joukowski profile cl error 
% minimization problem.
% 
% INPUTS: 
% - params(1), float: eta coordinate of the circle origin
% - params(2), float: xi coordinate of the circle origin
% - params(3), float: a parameter of the complex circle
% - alpha, float: alpha values to compute cl (rad)
% - cl, float: objective airfoil cl values
%
% OUTPUTS:
% - err, float: error between Joukowski cl curve and objective one
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 05/11/24, Boscariol Jacopo
%               Changes: release.

    % Joukowski cl computation
    cl_j = 2*pi*(params(2)/params(3) + alpha);

    % error computation
    % err = sum((cl - cl_j).^2);
    err = mean(abs(cl - cl_j));

end
