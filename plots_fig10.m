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

% adding vel magnitude contours
contourf(real(domain), imag(domain), vel_mag_c(:, :, idx), 30, ...
    'LineColor', 'none')
colorbar

% adding streamlines
h_stream = streamslice(real(domain), imag(domain), ...
    U_c(:, :, idx), V_c(:, :, idx), 3);
set(h_stream, 'Color', 'k')

title('Velocity Magnitude Around Joukowski Circle')
xlabel('Re(z)')
ylabel('Im(z)')

% velocity direction
figure
plot(real(zeta_circle), imag(zeta_circle), 'Color', 'k')
axis equal
hold on

% adding vel direction contours
contourf(real(domain), imag(domain), angles_c(:, :, idx), 30, ...
    'LineColor', 'none')
colorbar

% adding streamlines
h_stream = streamslice(real(domain), imag(domain), ...
    U_c(:, :, idx), V_c(:, :, idx), 3);
set(h_stream, 'Color', 'k')

title('Velocity Direction Around Joukowski Circle')
xlabel('Re(z)')
ylabel('Im(z)')

%% Joukowski profile plots

% velocity magnitude
figure
contourf(real(z_domain), imag(z_domain), vel_mag_p(:, :, idx), 30, ...
    'LineColor', 'none')
colorbar
hold on

% adding actual wing profile
plot(xW, yW, '-o', 'LineWidth', 1.5, 'MarkerSize', 1.5, ...
    'MarkerFaceColor', 'black', 'Color', 'black')

% adding actual tail profile
plot(xT, yT, '-o', 'LineWidth', 1.5, 'MarkerSize', 1.5, ...
    'MarkerFaceColor', 'black', 'Color', 'black')

title('Velocity Magnitude Around Joukowski Profile')
xlabel('x (m)')
ylabel('y (m)')

% applying suitable limits wrt the considered domain
xlim([-cT, L + 2*cT])
ylim([-cT, H + 2*cT])

% velocity direction
figure
contourf(real(z_domain), imag(z_domain), angles_p(:, :, idx), 30, ...
    'LineColor', 'none')
colorbar
hold on

% adding actual wing profile
plot(xW, yW, '-o', 'LineWidth', 1.5, 'MarkerSize', 1.5, ...
    'MarkerFaceColor', 'black', 'Color', 'black')

% adding actual tail profile
plot(xT, yT, '-o', 'LineWidth', 1.5, 'MarkerSize', 1.5, ...
    'MarkerFaceColor', 'black', 'Color', 'black')

title('Velocity Direction Around Joukowski Profile')
xlabel('x (m)')
ylabel('y (m)')

% applying suitable limits wrt the considered domain
xlim([-cT, L + 2*cT])
ylim([-cT, H + 2*cT])

%% Fig 10 plot

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
xlabel('\alpha (deg)')
ylabel('$\frac{c_L - c_{L,S}}{c_{L,S}}$ (-)', 'Interpreter', 'latex')
