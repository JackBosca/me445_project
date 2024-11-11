% Parameters NACA 5 digits LPSTT = 24012 (r, k1 obtained from tables)
r = 0.2025;  
k1 = 15.957;  
x = linspace(0, 1, 1000); 
yc = zeros(size(x));

U = 15;
c = 1;
L = 3*c;
H = c;

% Camberline definition For 0 < x < r
mask1 = (x > 0) & (x < r); 
yc(mask1) = (k1 / 6) * (x(mask1).^3 - 3*r*x(mask1).^2 + r^2 * (3 - r) * x(mask1));

% Case 2: For r < x < 1, apply the second equation
mask2 = (x >= r) & (x <= 1);
yc(mask2) = (k1 * r^3 / 6) * (1 - x(mask2));

% Camberline plot
figure;
plot(x, yc, 'LineWidth', 2);
xlabel('x',FontSize=14);
ylabel('y_c(x)',FontSize=14);
ylim([0 0.1]);
title('Camberline y_c(x)')

% Compute camberline gradient
theta = linspace(0, pi, 88); 
xi = (1 - cos(theta)) / 2; 
dyc_dx = zeros(size(xi));

% Case 1: For 0 < xi < r 
mask1 = (xi > 0) & (xi < r); 
dyc_dx(mask1) = (k1 / 6) * (3 * xi(mask1).^2 - 6 * r * xi(mask1) + r^2 * (3 - r));

% Case 2: For r < xi < 1 
mask2 = (xi >= r) & (xi <= 1);
dyc_dx(mask2) = -k1 * r^3 / 6;

% Fourier coefficients computation
A0 = 1/pi * trapz(theta, dyc_dx);
A1 = 2/pi * trapz(theta, dyc_dx .* cos(theta)); 

% Lift coefficients
alpha = linspace(-20, 30, 88) * pi / 180; 
Cl = 2 * pi * alpha + pi * (A1 - 2 * A0);

figure;
plot(alpha*180/pi, Cl);
xlabel('\alpha (deg)', FontSize=14);
ylabel('C_L',FontSize=14);
title('C_L VS \alpha')
% Downwash velocity 
Gamma = 0.5 * U * c * Cl;
d = sqrt(L^2 + H^2);
omega = Gamma / (2 * pi * d);
alpha_induced = atan2(omega, U);
U_eff = sqrt(U^2 + omega.^2);
alpha_eff = alpha - alpha_induced;

% Cl Tail wing
Cl_tail = 2 * pi * alpha_eff + pi * (A1 - 2 * A0);

% Relative CL for Fig 10
Relative_CL_tail = (Cl_tail - Cl) ./ Cl;

figure;
hold on;
plot(alpha * 180 / pi, Cl_tail);
plot(alpha * 180 / pi, Cl);
ylabel('C_L',FontSize=14)
xlabel('\alpha (deg)', FontSize=14)
title('Solo VS 2 airfoil configuration')
hold off;
legend('CL-tail', 'CL-wing');


% Convert the first and second columns to numeric arrays
alpha_experiment = load("data/fig10_tail_alpha.txt");
cl_experiment = load("data/fig10_tail_CL.txt");

figure;
hold on;
plot(alpha * 180 / pi, Relative_CL_tail);
scatter(alpha_experiment, cl_experiment);
ylabel('\Delta C_L',FontSize=16)
xlabel('\alpha (deg)', FontSize=16)
title('Fig 10 reproduction')
hold off;
legend('Numerical simulation', 'Paper data');

%Plot to check the computation of the CL
load("data/fig8_clW.mat");

figure;
hold on;
plot(fig8_clW(:, 1), fig8_clW(:, 2), 'Marker', 'p');
plot(alpha * 180 / pi, Cl);
hold off;
xlabel('\alpha (deg)',FontSize=15);
ylabel('CL_W',FontSize=15);


%Error analysis : residuals 
residuals = fig8_clW(:, 2)' - Cl;

figure;
hold on;

% Plot each residual line
for i = 1:length(fig8_clW(:, 1))
    % x-coordinate for the current residual point
    x_val = fig8_clW(i, 1); 
    
    % Residual value (y-coordinate of the end point of the line)
    y_val = residuals(i);
    
    % Plot a line from x-axis to the residual point
    plot([x_val, x_val], [0, y_val], 'k-', 'LineWidth', 1); % line from x-axis
    
    % Plot a dot at the top of each line
    plot(x_val, y_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 2.5); % dot at the end
end

% Labeling the plot
xlabel('\alpha (deg)',FontSize=15);
ylabel('Residuals (C_{l, experiment} - C_l)',FontSize=15);
title('Residuals of the C_L')
hold off;
