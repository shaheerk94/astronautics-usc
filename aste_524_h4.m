%% Imports

%% Constants 
R = 6378; % km, Earth Radius
u = 3.986e5; % km3/s2
g = 9.8/1000; % km/s^2
inverseB = 7.4; % km, scale height=B^-1
entry_interface = 122; % km, height of Earth atmosphere
e = exp(1); % euler's constant 
entry = entry_formulas();
orbit = orbital_mechanics_formulas();



%% Homework 4 Problems 

% Problem 1 
B = 0.1354; % km-1
ps = 1.225; % kg/m^3 
l_over_d = 1.05; % lift over drag
h = 25; % km, altitude 
Vf = 750/1000; %km/s
Ve = 7.5; %km/s
gamma_e = -1.2; % degrees 

% Part a 
s = entry.total_glided_range(R, l_over_d, Ve, g);
entry.disp_var('s', s)

% Part B 
a_g = entry.inst_glided_accel(Vf, g, R, l_over_d);
entry.disp_var('a_g', a_g)

% Part C 
theta = 45; % degrees
s_cross = entry.cross_ranging_dist_glided(R, l_over_d, theta); 
entry.disp_var('s_cross', s_cross) 

% Problem 2 
m = 2200; % kg
Cd = 1.5; 
l_over_d = 0.19; 
Ve = 7.5; % km/s

% Part A 
gamma_e = -5; % degrees 
amax = entry.max_accel_ballistic(B, Ve, gamma_e, e);
entry.disp_var('amax', amax);

B_m = 0.1354/1000; % m^-1 
S = 4.1; 

hcrit = entry.altitude_of_max_accel(B_m, S, Cd, m, ps, gamma_e);
entry.disp_var('hcrit', hcrit)

% Part B 
gamma_e = -2; % degrees 
Vf = 750/1000; % km/s 
a_g = entry.inst_glided_accel(Vf, g, R, l_over_d);
entry.disp_var('a_g', a_g);

% Part C 
g_km_s = 9.8/1000; % km/s 
s = entry.total_glided_range(R, l_over_d, Ve, g_km_s);
entry.disp_var('s', s);

% Part D 
theta = 45; % degrees
s_cross = entry.cross_ranging_dist_glided(R, l_over_d, theta);
entry.disp_var('s_cross', s_cross);

%% Problem 3 

% Part a 
amax = 8 * g; % km/s^2
Ve = 11; % km/s
e = exp(1);
B = 0.1354; % km-1
gamma_e = asin((2 * e * amax)/(-B * Ve^2)); 
gamma_e_deg = gamma_e * 180/pi;
entry.disp_var('gamma_e_deg', gamma_e_deg); 

