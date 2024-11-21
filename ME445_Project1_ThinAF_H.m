clear all; close all; clc;

% ---------- USER DEFINED ---------- 
V_inf = 25; 
alpha_deg = 5;
mu_x = -0.15; % Joukowski circle center (mu_x,mu_y)
mu_y = 0.096; 
a1 = [1.5,3,6]; % Tail location : L = 3*c = a1*c 
a2 = [0.5,1,2]; % Tail location : H = 1*c = a2*c 

%  ---------- EXP DATA ---------- 
alpha_exp = load('data/fig10_tail_alpha.txt');
Cl_exp = load('data/fig10_tail_CL.txt');
for i=1:length(alpha_exp)-1
    if alpha_exp(i) == -1.029185868000000
        alpha_exp(i)=[];
        Cl_exp(i)=[];
    end
    if alpha_exp(i) == 10.030721970000000
        alpha_slope2 = alpha_exp(i);
        Cl_slope2 = Cl_exp(i);
    end
    if alpha_exp(i) == -3.026113671000000
        alpha_slope1 = alpha_exp(i);
        Cl_slope1 = Cl_exp(i);
    end
end
slope_exp = (Cl_slope2-Cl_slope1)/(alpha_slope2-alpha_slope1);
b = Cl_slope2 - slope_exp*alpha_slope2;
f_slope_exp = slope_exp .* linspace(-5,10,30) + b;

% ---------- CIRCLE ----------  
alpha = deg2rad(alpha_deg); 
mu = mu_x + 1i * mu_y; 
R = sqrt((1 - mu_x)^2 + mu_y^2);
Gamma = 4 * pi * V_inf * R * sin(alpha + asin(mu_y / R)); % KUTTA CONDITION
% Gamma = 0;
theta = linspace(0, 2 * pi, 1000);
zeta_circle = R * exp(1i * theta) + mu;
chi_circle = real(zeta_circle);
eta_circle = imag(zeta_circle);
z_profil =  zeta_circle + 1 ./ zeta_circle; % JOUKOWSKI
x_profil = real(z_profil);
y_profil = imag(z_profil);
c = max(x_profil)-min(x_profil);

%  ---------- VELOCITY FIELD ---------- 
lim = 25;
subdivision = 450;
[xi, eta] = meshgrid(linspace(-lim, lim, subdivision), linspace(-lim, lim, subdivision));
zeta = xi + 1i * eta;
z = zeta + R^2 ./ zeta; 
dz_dzeta = 1 - (R^2 ./ zeta.^2); 
W_tilde = V_inf * exp(-1i * alpha) ...
    + (1i * Gamma) ./ (2 * pi * (zeta - mu)) ...
    - V_inf * R^2 * exp(1i * alpha) ./ (zeta - mu).^2;
W = W_tilde ./ dz_dzeta;

%  ---------- CIRCLE-VELOCITY
u_circle = real(W_tilde);
v_circle = -imag(W_tilde);
inside_circle = abs(zeta - mu) < R;
u_circle(inside_circle) = NaN;
v_circle(inside_circle) = NaN;

%  ---------- PROFIL-VELOCITY
u_profil = real(W);
v_profil = -imag(W);
u_profil(inside_circle) = NaN;
v_profil(inside_circle) = NaN;

% CIRCLE-FLOW
figure; hold on; grid on;
streamslice(xi, eta, [], u_circle, v_circle, [], [], [], [], 10);
plot(mu_x,mu_y,'r.','MarkerSize',10);
text(mu_x,mu_y-2,sprintf('(%.2f,%.2f)',mu_x,mu_y),...
'Color','Red','FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
plot(chi_circle, eta_circle, 'k-', 'LineWidth', 2); 
axis equal;
xlim([-5, 13]);
ylim([-5, 6]);
xlabel('$\chi$','Interpreter','latex');
ylabel('$\eta$','Interpreter','latex');

%  ---------- ANGLE AND AMPLITUDE @ TAIL LOCATION ---------- 
figure; hold on; grid on;
p1 = streamslice(real(z), imag(z), [], u_profil, v_profil, [], [],[],[],9); % Streamline 
p2 = plot(x_profil, y_profil, 'k-', 'LineWidth', 2); % Wing
for A1 = a1
    index_a1 = find(A1==a1);
    for A2 = a2
        index_a2 = find(A2==a2);
        L = A1.*c; H = A2.*c;
        x_tail = min(x_profil) + L.*cos(alpha);
        y_tail = c / 2 + L*sin(alpha);
        [~, idx_x] = min(abs(xi(1,:) - x_tail)); 
        [~, idx_y] = min(abs(eta(:,1) - y_tail)); 
        u_at_point = u_profil(idx_y, idx_x);
        v_at_point = v_profil(idx_y, idx_x);
        theta = atan2(v_at_point, u_at_point); 
        theta_deg = rad2deg(theta);
        amplitude = sqrt(u_at_point.^2 + v_at_point.^2);
        sprintf("a1 = %2f", "a2 = %.2f, x_T = %.2f, y_T = %.2f, θ = %.2f°, U = %.2f",A1,A2,x_tail,y_tail,theta_deg,amplitude);
        
        scale_factor = 0.5; 
        x_profil_scaled = x_tail + scale_factor .* x_profil; 
        y_profil_scaled = y_tail + scale_factor .* y_profil; 
        x_profil_scaled = x_profil_scaled + (scale_factor * c)/2; 
        
        % AIRFOIL-FLOW
        if A1 == 3 && A2 == 1
            theta_deg1 = theta_deg;
            p3 = plot(x_tail,y_tail,'r.', 'MarkerSize',10); % Location of tail
            plot(x_profil_scaled, y_profil_scaled, 'k', 'LineWidth',1.2); % Tail
            plot(x_profil, ones(length(x_profil)).*(min(y_profil)-0.5),'b','LineWidth',2); % Côte : Chord
            text(x_tail, y_tail+0.5, ...
                sprintf('θ = %.2f°\nU= %.2f', theta_deg, amplitude), ...
                'Color','red','FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); 
            text(x_tail,y_tail-0.50,sprintf('(%.2f,%.2f)',x_tail,y_tail),'FontSize',8,...
                'Color','red','HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            text(min(x_profil)+c/2, min(y_profil)-2.5, ...
                sprintf('c = %.2f\n α = %.1f°', c, alpha_deg), ...
                'Color','blue','FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); 
            axis equal;
            xlim([-5, 25]);
            ylim([-5, 7.5]);
            xlabel('$x$','Interpreter','latex');
            ylabel('$y$','Interpreter','latex');
        end
        if index_a1 == index_a2
            p3 = plot(x_tail,y_tail,'r.', 'MarkerSize',10); % Location of tail
            plot(x_profil_scaled, y_profil_scaled, 'Color',[0 0 0 0.3], 'LineWidth',1.2); % Tail
            text(x_tail, y_tail+0.5, ...
                sprintf('θ = %.2f°\nU= %.2f', theta_deg, amplitude), ...
                'Color','red','FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); 
            text(x_tail,y_tail-0.50,sprintf('(%.2f,%.2f)',x_tail,y_tail),'FontSize',8,...
                'Color','red','HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        end
        

        % ---------- TA Theory ---------- 
        theta_rad = deg2rad(theta_deg);
        
        % ---------- NACA23012 data ---------- 
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
        
        %  ---------- CAMBERLINE  ---------- 
        for i = 1:length(y_data)
            if y_data(i) > 0
                y_data_up   = [y_data_up,   y_data(i)];
                x_data_up   = [x_data_up,   x_data(i)];
            else
                y_data_down = [y_data_down, y_data(i)];
                x_data_down = [x_data_down, x_data(i)];
            end
        end
        x_naca      = interp1(x_data_up,   x_data_up,    x_data_down, 'linear');
        y_naca_up   = interp1(x_data_up,   y_data_up,    x_naca     , 'linear');
        y_naca_down = interp1(x_data_down, y_data_down,  x_naca     , 'linear');
        x_naca(1) = []; x_naca = [0 x_naca];
        y_naca_up(1) = []; y_naca_up = [0 y_naca_up];
        y_naca_down(1) = []; y_naca_down = [0 y_naca_down];
        for i = 1:length(x_naca)
            y_camber = [y_camber, (y_naca_up(i)+y_naca_down(i))/2];
        end
        
        %  ---------- INCLINED PROFILE ---------- 
        r = [cos(theta_rad), -sin(theta_rad); sin(theta_rad), cos(theta_rad)];
        rotated_up = r * [x_naca; y_naca_up];
        x_rot_up = rotated_up(1,:);
        y_rot_up = rotated_up(2,:);
        rotated_down = r * [x_naca; y_naca_down];
        x_rot_down = rotated_down(1,:);
        y_rot_down = rotated_down(2,:);
        rotated_camber = r * [x_naca; y_camber];
        x_rot_camber = rotated_camber(1,:);
        y_rot_camber = rotated_camber(2,:);
        
        %  ---------- TA Theory ---------- 
        if index_a1 == index_a2 
            dyc_dx = diff(y_camber) ./ diff(x_naca);
            dyc_dx_inclined = diff(y_rot_camber) ./ diff(x_naca);
            theta = acos(1 - 2.* x_naca); 
            dtheta = linspace(0,2*pi,31);
            % xi = 0.5 .* (1 - cos(theta));
            A0 = 1/pi .* trapz(dyc_dx,dtheta(:,2:end));
            for n = 1:10
                An = [An, 2/pi .* trapz(dyc_dx .* cos(n .* theta(:,2:end)), dtheta(2:end))];
            end
            A0_inclined = 1/pi .* trapz(dyc_dx_inclined,dtheta(:,2:end));
            for n = 1:10
                An_inclined = [An_inclined, 2/pi .* trapz(dyc_dx_inclined .* cos(n .* theta(:,2:end)), dtheta(2:end))];
            end
            alpha_TAF = linspace(-20,30,30);
            Cl_TAF = (2*pi.*alpha_TAF + pi*(An(1)-2*A0));
            Cl_TAF_inclined = (2*pi.*alpha_TAF + pi*(An_inclined(1)-2*A0_inclined));
            Cl_TAF_plot = (Cl_TAF - Cl_TAF_inclined) ./ Cl_TAF_inclined;

            
        end
    end
end
theta_deg = theta_deg1;


%  ---------- ERROR ---------- 
% error = abs(Cl_TAF_plot-f_slope_exp);

%  ---------- FIGURES ---------- 


% AIRFOIL
figure; hold on; grid on;
p1 = plot(x_naca, y_naca_up,  'b-', 'LineWidth', 1, 'Color',[0 0 0 0.2]);
plot(x_naca, y_naca_down, 'b-', 'LineWidth', 1, 'Color',[0 0 0 0.2]); 
p2 = plot(x_naca, y_camber, 'r-', 'LineWidth', 1, 'Color',[1 0 0 0.2]);
p3 = plot(x_rot_up, y_rot_up, 'k-', 'LineWidth', 1); % Rotated upper surface
plot(x_rot_down, y_rot_down, 'k-', 'LineWidth', 1); % Rotated lower surface
p4 = plot(x_rot_camber, y_rot_camber, 'r-', 'LineWidth', 1); % Rotated camber line
p5 = plot(x_naca(1,3:end), dyc_dx(1,2:end), 'b-', 'Color',[0 0 1 0.4]);
p6 = plot(x_naca(1,3:end), dyc_dx_inclined(1,2:end), 'b-');
axis equal;
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
legend([p3, p4, p6], 'NACA23012', '$y_{c}$', ...
    '$dy_c/dx$', 'Interpreter', 'latex');

% (CL-CLS)/CLS
figure; hold on; grid on;
p1 = plot(alpha_TAF,Cl_TAF_plot,'b');
p2 = plot(alpha_exp,Cl_exp,'ko','MarkerSize',5);
plot(-1.029185868000000,0.3186,'ko');
p3 = plot(linspace(-5,10,30),f_slope_exp,'Color',[1 0 0 0.4]);
xlabel('$\alpha$','Interpreter','latex');
ylabel('$(C_L - C_{L,S})/C_{L,S}$','Interpreter','latex');
legend([p1,p2,p3], 'TA Theory','Experimental','Slope approx.','Interpreter','Latex');
