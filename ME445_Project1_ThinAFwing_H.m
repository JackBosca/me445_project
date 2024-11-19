close all; clear all; clc;

% -------------------- User-defined angle of attack
alpha_deg = 0;
alpha_rad = deg2rad(alpha_deg);

% -------------------- Load NACA data
data = load('data/NACA23012.txt');
x_data = data(:,1); 
y_data = data(:,2);
x_data_up = [];
x_data_down = [];
y_data_up = [];
y_data_down = [];
y_camber = [];
An = [];
An_inclined = [];

% -------------------- Separation of y_up and y_down
for i = 1:length(y_data)
    if y_data(i) > 0
        y_data_up   = [y_data_up,   y_data(i)];
        x_data_up   = [x_data_up,   x_data(i)];
    else
        y_data_down = [y_data_down, y_data(i)];
        x_data_down = [x_data_down, x_data(i)];
    end
end

% -------------------- Setting up y_up and y_down for same x values
x_naca      = interp1(x_data_up,   x_data_up,    x_data_down, 'linear');
y_naca_up   = interp1(x_data_up,   y_data_up,    x_naca     , 'linear');
y_naca_down = interp1(x_data_down, y_data_down,  x_naca     , 'linear');
x_naca(1) = []; x_naca = [0 x_naca];
y_naca_up(1) = []; y_naca_up = [0 y_naca_up];
y_naca_down(1) = []; y_naca_down = [0 y_naca_down];

% -------------------- Camber line
for i = 1:length(x_naca)
    y_camber = [y_camber, (y_naca_up(i)+y_naca_down(i))/2];
end

% -------------------- Apply rotation for angle of attack
R = [cos(alpha_rad), -sin(alpha_rad); sin(alpha_rad), cos(alpha_rad)];

% Upper surface
rotated_up = R * [x_naca; y_naca_up];
x_rot_up = rotated_up(1,:);
y_rot_up = rotated_up(2,:);

% Lower surface
rotated_down = R * [x_naca; y_naca_down];
x_rot_down = rotated_down(1,:);
y_rot_down = rotated_down(2,:);

% Camber line
rotated_camber = R * [x_naca; y_camber];
x_rot_camber = rotated_camber(1,:);
y_rot_camber = rotated_camber(2,:);

% -------------------- Thin airfoil theory
dyc_dx = diff(y_camber) ./ diff(x_naca);
dyc_dx_inclined = diff(y_rot_camber) ./ diff(x_naca);

theta = acos(1 - 2.* x_naca); 
dtheta = linspace(0,2*pi,31);
xi = 0.5 .* (1 - cos(theta));

A0 = 1/pi .* trapz(dyc_dx,dtheta(:,2:end));
for n = 1:10
    An = [An, 2/pi .* trapz(dyc_dx .* cos(n .* theta(:,2:end)), dtheta(2:end))];
end

A0_inclined = 1/pi .* trapz(dyc_dx_inclined,dtheta(:,2:end));
for n = 1:10
    An_inclined = [An_inclined, 2/pi .* trapz(dyc_dx_inclined .* cos(n .* theta(:,2:end)), dtheta(2:end))];
end

alpha_TAF = linspace(-20,30,30);
Cl_TAF = (2*pi.*alpha_TAF + pi*(An(1)-2*A0))/100;
Cl_TAF_inclined = (2*pi.*alpha_TAF + pi*(An_inclined(1)-2*A0_inclined))/100;

Cl_TAF_plot = (Cl_TAF - Cl_TAF_inclined) ./ Cl_TAF_inclined;

% -------------------- Exp data
alpha_exp = load('data/fig10_wing_alpha.txt');
Cl_exp = load('data/fig10_wing_CL.txt');
for i=1:length(alpha_exp)
    if alpha_exp(i) == -1.029185868000000
        alpha_exp(i)=[];
        Cl_exp(i)=[];
    end
    if alpha_exp(i) == 30
        alpha_slope2 = alpha_exp(i);
        Cl_slope2 = Cl_exp(i);
    end
    if alpha_exp(i) == -18.003072200000000
        alpha_slope1 = alpha_exp(i);
        Cl_slope1 = Cl_exp(i);
    end
end

slope_exp = (Cl_slope2-Cl_slope1)/(alpha_slope2-alpha_slope1);
b = Cl_slope2 - slope_exp*alpha_slope2;
f_slope_exp = slope_exp .* linspace(-20,30,30) + b;

% -------------------- Error
error = abs(Cl_TAF_plot-f_slope_exp);

% -------------------- Figure Plot
figure; hold on; grid on;
p1 = plot(x_naca, y_naca_up,  'b--', 'LineWidth', 1, 'Color',[0 0 0 0.2]);
plot(x_naca, y_naca_down, 'b--', 'LineWidth', 1, 'Color',[0 0 0 0.2]); 
p2 = plot(x_naca, y_camber, 'r--', 'LineWidth', 1, 'Color',[1 0 0 0.2]);
p3 = plot(x_rot_up, y_rot_up, 'k-', 'LineWidth', 1); % Rotated upper surface
plot(x_rot_down, y_rot_down, 'k-', 'LineWidth', 1); % Rotated lower surface
p4 = plot(x_rot_camber, y_rot_camber, 'r-', 'LineWidth', 1); % Rotated camber line
p5 = plot(x_naca(1,3:end), dyc_dx(1,2:end), 'b--', 'Color',[0 0 1 0.4]);
p6 = plot(x_naca(1,3:end), dyc_dx_inclined(1,2:end), 'b-');
axis equal;
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
legend([p3, p4, p6], 'NACA23012', '$y_{c}$', ...
    '$dy_c/dx$', 'Interpreter', 'latex');

figure; hold on; grid on;
p1 = plot(alpha_TAF,Cl_TAF_plot,'b');
p2 = plot(alpha_exp,Cl_exp,'ko','MarkerSize',5);
p3 = plot(linspace(-20,30,30),f_slope_exp,'Color',[1 0 0 0.4]);
xlabel('$\alpha$','Interpreter','latex');
ylabel('$(C_L - C_{L,S})/C_{L,S}$','Interpreter','latex');
legend([p1,p2,p3], 'TA Theory','Experimental','Slope approx.','Interpreter','Latex');
