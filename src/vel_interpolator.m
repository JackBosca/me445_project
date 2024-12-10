function [U_LH, V_LH, alpha_eff] = vel_interpolator(z_domain, U, V, L, H)
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
% - U_LH, float: x-velocities at query point for each AoA
% - V_LH, float: y-velocities at query point for each AoA
% - alpha_eff, float: effective alpha at query point for each AoA
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

    % computing step sizes along rows (x) and columns (y)
    dx = diff(real(z_domain), 1, 2);
    dy = diff(imag(z_domain), 1, 1);

    % finding mean step sizes in both directions
    mean_dx = mean(abs(dx(:)));
    mean_dy = mean(abs(dy(:)));

    % defining the search radius as twice the largest mean step size
    search_radius = 2*max(mean_dx, mean_dy);

    % finding indices of points near (L, H)
    distances = sqrt((x - L).^2 + (y - H).^2);
    nearby_idx = distances <= search_radius;

    % defining close points to optimize interpolation time
    x_near = x(nearby_idx);
    y_near = y(nearby_idx);

    % output init
    U_LH = zeros(n_alpha, 1);
    V_LH = zeros(n_alpha, 1);

    % main loop
    for ii = 1:n_alpha
        % selecting current u, v near (L, H)
        u = U(:, :, ii);
        v = V(:, :, ii);
        u = u(nearby_idx);
        v = v(nearby_idx);

        % interpolating using scatteredInterpolant()
        F_U = scatteredInterpolant(x_near, y_near, u, 'linear');
        U_LH(ii) = F_U(L, H);

        F_V = scatteredInterpolant(x_near, y_near, v, 'linear');
        V_LH(ii) = F_V(L, H);
    end

    % computing angles at (L, H)
    alpha_eff = atan2(V_LH, U_LH);

end
