function [U_LH, V_LH, angles] = vel_interpolator(z_domain, U, V, L, H)
% vel_interpolator - interpolates the velocity field (U, V) in the point
% (L, H) of the considered geometrical domani z_domain (which does not need
% to be generated through meshgrid).
% 
% INPUTS: 
% - z_domain, complex float: domain in the Joukowski profile plane
% - U, float: x-velocity field for each angles of attack
% - V, float: y-velocity field for each angles of attack
% - L, float: x-coordinate of the query point
% - H, float: y-coordinate of the query point
%
% OUTPUTS:
% - U_LH, float: x-velocities at query point for each angles of attack
% - V_LH, float: y-velocities at query point for each angles of attack
% - angles, float: velocity angles at query point for each angles of attack
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 26/11/24, Boscariol Jacopo
%               Changes: release.
    
    % extracting sizes
    n_alpha = size(U, 3);

    % flattening z_domain
    x = real(z_domain(:));
    y = imag(z_domain(:));
    
    % output init
    U_LH = zeros(n_alpha, 1);
    V_LH = zeros(n_alpha, 1);

    % main loop
    for ii = 1:n_alpha
        % selecting current u, v
        u = U(:, :, ii);
        v = V(:, :, ii);

        % interpolating using scatteredInterpolant()
        F_U = scatteredInterpolant(x, y, u(:), 'linear');
        U_LH(ii) = F_U(L, H);

        F_V = scatteredInterpolant(x, y, v(:), 'linear');
        V_LH(ii) = F_V(L, H);
    end

    % computing angles at (L, H)
    angles = atan2(V_LH, U_LH);

end
