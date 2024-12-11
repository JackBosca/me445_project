clearvars -except Cl_TAT_rot cls
close all; clc;

% ---------- USER DEFINED ---------- 
V_inf = 25; 
alpha_deg = 0.4;
alpha_rad = deg2rad(alpha_deg);
a1 = [3]; % Tail location : L = 3*c = a1*c 
a2 = [1]; % Tail location : H = 1*c = a2*c 

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


% ---------------------- DEF-CIRCLE ----------------------
mu_y = 0.085; 
mu_x = 0.025; 
mu   = mu_y + 1i*mu_x;
R    = abs(mu - 1); 
theta_rad = linspace(0,2*pi,1000);
zeta_circle = R * exp(1i * theta_rad) + mu;
chi_circle = real(zeta_circle);
eta_circle = imag(zeta_circle);

% ---------------------- DEF-PROFILE ----------------------
z_profil = zeta_circle + 1 ./ zeta_circle;
x_profil = real(z_profil); y_profil = imag(z_profil);
c_profil = max(x_profil)-min(x_profil);
scaling_factor = 50/c_profil;

% Rescale of velocity : 
V_inf = V_inf/scaling_factor;

% ---------------------- VELOCITY FIELD ----------------------
Gamma = 4 * pi * V_inf * R * sin(alpha_rad + asin(mu_y / R)); % KUTTA CONDITION

lim = 25;
subdivision = 450;
[xi, eta] = meshgrid(linspace(-lim, lim, subdivision), linspace(-lim, lim, subdivision));
zeta = xi + 1i * eta;
z = zeta + R^2 ./ zeta; 
dz_dzeta = 1 - (R^2 ./ zeta.^2); 
W_tilde = V_inf * exp(-1i * alpha_rad) ...
    + (1i * Gamma) ./ (2 * pi * (zeta - mu)) ...
    - V_inf * R^2 * exp(1i * alpha_rad) ./ (zeta - mu).^2;
W = W_tilde ./ dz_dzeta;

u_circle =  real(W_tilde);
v_circle = -imag(W_tilde);
inside_circle = abs(zeta - mu) < R;

u_profil = real(W);
v_profil = -imag(W);

% ---------------------- TAIL LOCATION ----------------------
L = a1.*c_profil; H = a2.*c_profil;
x_tail_pos = min(x_profil) + L.*cos(alpha_rad);
y_tail_pos = c_profil / 2 + L*sin(alpha_rad);
[~, idx_x] = min(abs(xi(1,:) - x_tail_pos)); 
[~, idx_y] = min(abs(eta(:,1) - y_tail_pos)); 
u_at_point = u_profil(idx_y, idx_x);
v_at_point = v_profil(idx_y, idx_x);

theta_rad = atan2(v_at_point, u_at_point);
theta_deg = rad2deg(theta_rad);
amplitude = sqrt(u_at_point.^2 + v_at_point.^2);

% ---------------------- TA THEORY ----------------------
% ---------------------- Camberline (Joukowski Profile)
x_profil_up = [];
x_profil_down = [];
y_profil_up = [];
y_profil_down = [];
y_camber = [];
for i = 1:length(y_profil)
    if y_profil(i) > 0
        y_profil_up   = [y_profil_up,   y_profil(i)];
        x_profil_up   = [x_profil_up,   x_profil(i)];
    else
        y_profil_down = [y_profil_down, y_profil(i)];
        x_profil_down = [x_profil_down, x_profil(i)];
    end
end

x_profil_interp      = interp1(x_profil_up,   x_profil_up,    x_profil_down, 'linear');
y_profil_up   = interp1(x_profil_up,   y_profil_up,    x_profil_interp     , 'linear');
y_profil_down = interp1(x_profil_down, y_profil_down,  x_profil_interp     , 'linear');

for i = 1:length(x_profil_interp)
    y_camber = [y_camber, (y_profil_up(i)+y_profil_down(i))/2];
end

% ---------------------- Camberline Rotation
r = [cos(theta_rad), -sin(theta_rad); sin(theta_rad), cos(theta_rad)];
rotated_up = r * [x_profil_interp; y_profil_up];
x_rot_up = rotated_up(1,:);
y_rot_up = rotated_up(2,:);
rotated_down = r * [x_profil_interp; y_profil_down];
x_rot_down = rotated_down(1,:);
y_rot_down = rotated_down(2,:);
rotated_camber = r * [x_profil_interp; y_camber];
x_rot_camber = rotated_camber(1,:);
y_rot_camber = rotated_camber(2,:);

% Rescale of tail (Paper : c_T = c_W/2)
x_profil_interp = x_profil_interp./2;
y_profil_up = y_profil_up./2;
y_profil_down = y_profil_down./2;

x_rot_up = x_rot_up./2;
y_rot_up = y_rot_up./2;
x_rot_down = x_rot_down./2;
y_rot_down = y_rot_down./2;

x_rot_camber = x_rot_camber./2;
y_rot_camber = y_rot_camber./2;

% ---------------------- TA Theory formulas
dyc_dx = diff(y_camber)./diff(x_profil_interp);
dyc_dx_rot = diff(y_rot_camber)./diff(x_rot_camber);

theta = acos(1-2.*x_profil_interp);
theta_rot = acos(1-2.*x_rot_camber);
dtheta = linspace(0,2*pi,407);

A0     = 1/pi .* trapz(dyc_dx     ,dtheta);
A0_rot = 1/pi .* trapz(dyc_dx_rot ,dtheta);
A1     = 2/pi .* trapz(dyc_dx     .*cos(theta(:,2:end))       , dtheta);
A1_rot = 2/pi .* trapz(dyc_dx_rot .*cos(theta_rot(:,2:end))   , dtheta); 
A2     = 2/pi .* trapz(dyc_dx     .*cos(2.*theta(:,2:end))    , dtheta);
A2_rot = 2/pi .* trapz(dyc_dx_rot .*cos(2.*theta_rot(:,2:end)), dtheta);

alpha_TAT = linspace(-20,30,30);
Cl_TAT     = 2*pi.*alpha_TAT + pi*(A1-2*A0);
Cl_TAT_rot = 2*pi.*alpha_TAT + pi*(A1_rot-2*A0_rot);
Cl_TAT = (Cl_TAT - Cl_TAT_rot) ./ Cl_TAT_rot;

CL = 2*pi*theta_rad + pi*(A1_rot - 2*A0_rot);

% Results
figure; hold on; grid on;
plot(alpha_TAT,Cl_TAT,'b')
plot(alpha_exp,Cl_exp,'ko')
axis square;

