close all; clear all; clc;

data = load('NACA23012.txt');
x_data = data(:,1); 
y_data = data(:,2);
x_data_up = [];
x_data_down = [];
y_data_up = [];
y_data_down = [];
y_camber = [];
An = [];

% --------------------  Separation of y_up and y_down
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

% -------------------- Thin airfoil theory
dyc_dx = diff(y_camber) ./ diff(x_naca);
theta = acos(1 - 2.* x_naca); 
xi = 0.5 .* (1 - cos(theta));
A0 = 1/pi .* trapz(dyc_dx,theta(:,2:end));
for n = 1:10
    An = [An, 2/pi .* trapz(dyc_dx .* cos(n .* theta(:,2:end)), theta(2:end))]
end

alpha = linspace(-5,10,30)
Cl = (2*pi.*alpha + pi*(An(1)-2*A0))/100

% -------------------- Exp data
alpha_exp = load('fig10_tail_alpha.txt')
Cl_exp = load('fig10_tail_CL.txt')
for i=1:length(alpha)
    if alpha_exp(i) == -1.029185868000000
        alpha_exp(i)=[];
        Cl_exp(i)=[];
    end
    if alpha_exp(i) == 10.030721970000000
        alpha_slope2 = alpha_exp(i);
        Cl_slope2 = Cl_exp(i);
    end
    if alpha_exp(i) == -3.026113671000000
        alpha_slope1 = alpha_exp(i)
        Cl_slope1 = Cl_exp(i)
    end
end

slope_exp = (Cl_slope2-Cl_slope1)/(alpha_slope2-alpha_slope1);
b = Cl_slope2 - slope_exp*alpha_slope2;
f_slope = slope_exp .* linspace(-5,10,30) + b

% -------------------- Figure Plot
figure; hold on; grid on;
p1 = plot(x_naca, y_naca_up,  'b.-', 'LineWidth', 1) ;
plot(x_naca,y_naca_down, 'b.-', 'LineWidth', 1); 
% plot(x_naca,y_naca_down, 'b.-', 'LineWidth', 1);
p2 = plot(x_naca,y_camber, 'r.-');
p3 = plot(x_naca(1,3:end), dyc_dx(1,2:end), '-');
axis equal;
xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');
legend([p1,p2,p3],'NACA 23012','Camber line','$dy_c/dx$','Interpreter','latex');

figure; hold on; grid on;
p1 = plot(alpha,Cl,'b')
p2 = plot(alpha_exp,Cl_exp,'ko-','MarkerSize',5)
plot(-1.029185868000000,0.3186,'ko')
plot([alpha_slope1 alpha_slope2],[Cl_slope1 Cl_slope2],'g')
plot(linspace(-5,10,30),f_slope)
xlabel('$\alpha$','Interpreter','latex');
ylabel('$(C_L - C_{L,S})/C_{L,S}$','Interpreter','latex');
legend([p1,p2], 'Thin airfoil theory','Experimental','Interpreter','Latex')
