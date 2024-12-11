%{
plots_fig10.m - Plots script for the Aerodynamics ME-445 project, figure 10. 
Selected paper: https://doi.org/10.1007/s00348-017-2429-4.

CALLED FUNCTIONS: -

CALLED DATA FILES: -

REVISIONS:
- #v0, 26-11-2024, Release, Boscariol Jacopo

Changes: -
%}

% set colors
col_geom = 'blue';
colW = 'red';
colT = 'green';

%% Joukowski circle plots

% velocity magnitude
figure
plot(real(zeta_circle), imag(zeta_circle), 'Color', 'k')
axis equal
hold on
plot(opt_params(1), opt_params(2), 'kx')

% adding vel magnitude contours
contourf(real(domain), imag(domain), vel_mag_c(:, :, idx), 30, ...
    'LineColor', 'none')
colormap('jet')
cb = colorbar;
cb.Label.String = 'Velocity Magnitude (m/s)';

% adding streamlines
h_stream = streamslice(real(domain), imag(domain), ...
    U_c(:, :, idx), V_c(:, :, idx), 3);
set(h_stream, 'Color', 'k')

title(sprintf(['Velocity Magnitude, Joukowski Circle - $\\alpha = %.2f^\\circ$, ' ...
    '$U_{\\infty} = %.2f$ m/s'], rad2deg(alphaT(idx)), U_inf), 'Interpreter', 'latex');
xlabel('Re(z)')
ylabel('Im(z)')

% velocity direction
figure
plot(real(zeta_circle), imag(zeta_circle), 'Color', 'k')
axis equal
hold on
plot(opt_params(1), opt_params(2), 'kx')

% adding vel direction contours
contourf(real(domain), imag(domain), angles_c(:, :, idx), 30, ...
    'LineColor', 'none')
colormap('jet')
cb = colorbar;
cb.Label.String = 'Velocity Direction (deg)';

% adding streamlines
h_stream = streamslice(real(domain), imag(domain), ...
    U_c(:, :, idx), V_c(:, :, idx), 3);
set(h_stream, 'Color', 'k')

title(sprintf(['Velocity Direction, Joukowski Circle - $\\alpha = %.2f^\\circ$, ' ...
    '$U_{\\infty} = %.2f$ m/s'], rad2deg(alphaT(idx)), U_inf), 'Interpreter', 'latex');
xlabel('Re(z)')
ylabel('Im(z)')

%% Joukowski profile plots

% velocity magnitude
figure
contourf(real(z_domain), imag(z_domain), vel_mag_p(:, :, idx), 30, ...
    'LineColor', 'none')
colormap('jet')
cb = colorbar;
cb.Label.String = 'Velocity Magnitude (m/s)';
hold on

% adding actual wing profile
plot(xW, yW, '-', 'LineWidth', 1.5, 'MarkerSize', 1.5, ...
    'MarkerFaceColor', 'black', 'Color', 'black')

% adding actual tail profile
plot(xT, yT, '-', 'LineWidth', 1.5, 'MarkerSize', 1.5, ...
    'MarkerFaceColor', 'black', 'Color', 'black')

% adding text annotation at tail location
vel_mag = norm([U_LH(idx), V_LH(idx)]);
text_string = sprintf('Velocity magnitude = %.2f m/s,\n Velocity direction = %.2f°', ...
    vel_mag, rad2deg(alpha_effT(idx)));
text(L, H, text_string, 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right');

title(sprintf(['Velocity Magnitude, Joukowski Profile - $\\alpha = %.2f^\\circ$, ' ...
    '$U_{\\infty} = %.2f$ m/s'], rad2deg(alphaT(idx)), U_inf), 'Interpreter', 'latex');
xlabel('x (m)')
ylabel('y (m)')

% applying suitable limits wrt the considered domain
xlim([-cT, L + 2*cT])
ylim([-cT, H + 2*cT])

% velocity direction
figure
contourf(real(z_domain), imag(z_domain), angles_p(:, :, idx), 30, ...
    'LineColor', 'none')
colormap('jet')
cb = colorbar;
cb.Label.String = 'Velocity Direction (deg)';
hold on

% adding actual wing profile
plot(xW, yW, '-', 'LineWidth', 1.5, 'MarkerSize', 1.5, ...
    'MarkerFaceColor', 'black', 'Color', 'black')

% adding actual tail profile
plot(xT, yT, '-', 'LineWidth', 1.5, 'MarkerSize', 1.5, ...
    'MarkerFaceColor', 'black', 'Color', 'black')

% adding text annotation at tail location
vel_mag = norm([U_LH(idx), V_LH(idx)]);
text_string = sprintf('Velocity magnitude = %.2f m/s,\n Velocity direction = %.2f°', ...
    vel_mag, rad2deg(alpha_effT(idx)));
text(L, H, text_string, 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right');

title(sprintf(['Velocity Direction, Joukowski Profile - $\\alpha = %.2f^\\circ$, ' ...
    '$U_{\\infty} = %.2f$ m/s'], rad2deg(alphaT(idx)), U_inf), 'Interpreter', 'latex');
xlabel('x (m)')
ylabel('y (m)')

% applying suitable limits wrt the considered domain
xlim([-cT, L + 2*cT])
ylim([-cT, H + 2*cT])

%% Fig 10 plot

% filtering ratio, alpha_effT not meaningful values (threshold = 1)
ratio = ratio(abs(ratio) < 1);
alpha_effT = alpha_effT(abs(ratio) < 1);

figure
% adding computed ratio with filtered alphaT
plot(rad2deg(alphaT(abs(ratio) < 1)), ratio, 'k^-')
hold on
grid on

% if L, H are default, plot the paper data
if L == 3*c && H == c
    plot(fig10_T(:, 1), fig10_T(:, 2), 'ko')
    legend('Thin Airfoil Theory', 'Paper data')
end

title('Fig. 10 - Lift Interaction Parameter')
xlabel('$\alpha$ (deg)', 'Interpreter', 'latex')
ylabel('$\frac{c_L - c_{L,S}}{c_{L,S}}$ (-)', 'Interpreter', 'latex')

%% Fig 9 plot

if L == 3*c && H == c
    % search linear interval in clT(alpha) curve
    idx1 = find(fig10_T(:, 1) > alphas_nT, 1);
    idx2 = find(fig10_T(:, 1) > alphas_pT, 1);
    alphaT_fig9 = alphaT(idx1:idx2);
    clT_fig9 = clT(idx1:idx2);

    figure
    plot(fig9_clT(:, 1), fig9_clT(:, 2), 'Marker', 'square', 'Color', 'black')
    hold on
    grid on
    xlabel('\alpha (deg)')
    ylabel('c_l (-)')

    % thin airfoil theory
    plot(rad2deg(alphaT_fig9), clT_fig9, 'Marker', 'o', 'MarkerSize', 3, ...
        'Color', 'red')

    % adding errorbars
    errorbar(rad2deg(alphaT_fig9), clT_fig9, ...
        uyT, uyT, rad2deg(uxT), rad2deg(uxT), 'CapSize', 3, 'Color', 'red')

    legend('Paper data', 'Thin Airfoil Theory')
    title('Fig. 9 - Lift Coefficient in 2-Airfoils Configuration')
end
