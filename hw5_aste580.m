%% ASTE 580 HW 5
% Shaheer Khan
% 4/9/2025
RE = 6400; %km
u = 4e5; % km^3/s^2
J2 = 1e-3;
function disp_var(var, value)
    disp([var, ' = ', num2str(value)]);
end

%% Problem 1
disp('-------------')
disp('Problem 1')


nE = 2 * pi/(365.25 * 24 * 3600); 
result = zeros(size(r));
disp_var('nE', nE);

am = (1.5 * sqrt(u) * RE^2 * J2/nE)^(2/7);
disp_var('am', am)

%% Part b
a = 6500:10:am; % km

i = acos(-(a./am).^(7/2));
i_deg = rad2deg(i);
plot(a, i_deg, 'b', 'LineWidth', 2)
xlabel('Orbit Radius a (km)')
ylabel('Inclination (degrees)')
title('Sun-Synchronous Orbit Inclination vs Orbit Radius')
grid on

% Part c 

a = (-am^(7/2) * cosd(97))^(2/7);

T = 2 * pi * sqrt(a^3/u);
disp_var('T @ i=97', T);

Tm = 2 * pi * sqrt(am^3/u);
disp_var('Tm @ i=180', Tm);

%% Problem 2 
% Part a 
% Constants
mu = 4e5; % km^3/s^2
RE = 6400; % km
J2 = 1e-3;
T = 43200; % seconds (12 hours)
i = asin(2/sqrt(5)); % radians
rp = 7000; % km

a = (mu * T^2 / (4 * pi^2))^(1/3);
e = 1 - rp / a;
p = a * (1 - e^2);
n = 2 * pi / T;
Omega_dot = -1.5 * J2 * (RE^2 / p^2) * n * cos(i); % rad/s
Omega_dot_deg_day = Omega_dot * (180/pi) * 86400;
t_1deg_days = 1 / abs(Omega_dot_deg_day);

disp_var('omega_dot_rad_s', Omega_dot);
disp_var('omega_dot_deg_day', Omega_dot_deg_day);
disp_var('Time to precess 1 degree (days)', t_1deg_days)

% Part b 
E = (e + cosd(90))/(1 + e*cosd(90));
tsouth = sqrt(a^3/u) * (E - e*sin(E));
tnorth = T - 2 * tsouth;
disp_var('tnorth', tnorth);

% Part c=
i = 63.455 - 0.5; 
wdot = 1.5 * n * RE^2 * J2/p^2 * (2 - 5/2 *(sind(i)^2));
disp_var('wdot', wdot);
t = pi/wdot; 
disp_var('t', t);
