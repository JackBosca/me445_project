clear all; close all; clc;

% CIRCLE
mu_x = 0.05; 
mu_y = 0.3;  
mu = mu_x + 1i * mu_y; 
V_inf = 30; 
alpha = deg2rad(5); 
R = sqrt((1 - mu_x)^2 + mu_y^2);
Gamma = 4 * pi * V_inf * R * sin(alpha + asin(mu_y / R)); % KUTTA CONDITION
% Gamma = 0;
theta = linspace(0, 2 * pi, 1000);
zeta_circle = R * exp(1i * theta) + mu;
chi_circle = real(zeta_circle);
eta_circle = imag(zeta_circle);

z_profil =  zeta_circle + R^2 ./ zeta_circle; % JOUKOWSKI
x_profil = real(z_profil);
y_profil = imag(z_profil);

c = max(x_profil)-min(x_profil);

% VELOCITY FIELD
lim = 12;
[xi, eta] = meshgrid(linspace(-lim, lim, 200), linspace(-lim, lim, 200));
zeta = xi + 1i * eta;
z = zeta + R^2 ./ zeta; 
dz_dzeta = 1 - (R^2 ./ zeta.^2); 
W_tilde = V_inf * exp(-1i * alpha) ...
    + (1i * Gamma) ./ (2 * pi * (zeta - mu)) ...
    - V_inf * R^2 * exp(1i * alpha) ./ (zeta - mu).^2;
W = W_tilde ./ dz_dzeta;

% CIRCLE-VELOCITY
u_cercle = real(W_tilde);
v_cercle = -imag(W_tilde);

% PROFIL-VELOCITY
u_profil = real(W);
v_profil = -imag(W);

% ANGLE AND AMPLITUDE @ TAIL LOCATION
x_target = min(x_profil) + 3 * c;
y_target = c / 2;
[~, idx_x] = min(abs(xi(1,:) - x_target)); % Trouver l'indice de x le plus proche
[~, idx_y] = min(abs(eta(:,1) - y_target)); % Trouver l'indice de y le plus proche
u_at_point = u_profil(idx_y, idx_x);
v_at_point = v_profil(idx_y, idx_x);
theta = atan2(v_at_point, u_at_point); % atan2 gère les signes et donne l'angle en radians
theta_deg = rad2deg(theta);
amplitude = sqrt(u_at_point^2 + v_at_point^2);

% FIGURES
% CIRCLE-PLAN
figure; hold on; grid on;
streamslice(xi, eta, u_cercle, v_cercle);
plot(mu_x,mu_y,'r.','MarkerSize',10);
text(mu_x,mu_y-R-0.75,sprintf('(%.2f,%.2f)',mu_x,mu_y),...
    'Color','Red','FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
plot(chi_circle, eta_circle, 'k-', 'LineWidth', 2); 
axis equal;
xlim([-5, lim]);
ylim([-5, lim]);
xlabel('$\chi$','Interpreter','latex');
ylabel('$\eta$','Interpreter','latex');

% AIRFOIL-PLAN
figure; hold on; grid on;
p1 = streamslice(real(z), imag(z), u_profil, v_profil);
p2 = plot(x_profil, y_profil, 'k-', 'LineWidth', 2);
p3 = plot(x_target,y_target,'k.', 'MarkerSize',8);
text(x_target, y_target+0.25, ...
    sprintf('θ = %.2f°\nU= %.2f', theta_deg, amplitude), ...
    'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); 
text(x_target,y_target-0.25,sprintf('(%.2f,%.2f)',x_target,y_target),'FontSize',8,...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
axis equal;
xlim([-5, lim]);
ylim([-5, lim]);
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
