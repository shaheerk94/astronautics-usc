%% Constants 
R = 6378; % km, Earth Radius
u = 3.986e5; % km3/s2
g_ms = 9.8; %m/s^2
g_kms = 9.8/1000; % km/s^2
entry_interface = 122; % km, height of Earth atmosphere

entry = entry_formulas();
orbit = orbital_mechanics_formulas();

function disp_var(var, value)
    disp([var, ' = ', num2str(value)]);
end

%% ASTE 524 midterm

%% Problem 2b 

thetai = 110 * pi/180; % radians
r_earth = 1.4952e8; % km
r_uranus = 2.8696e9; % km 

P_earth = 365 * 86400; % s
P_uranus = 30687 * 86400; % s


usun = 1.327e11; %km^3/s^2 

% earth = chaser 
wc = orbit.angular_speed(usun, r_earth); 
disp_var('wc', wc);
wt = orbit.angular_speed(usun, r_uranus); 
disp_var('wt', wt); 
thetaf = orbit.phase_angle_rendezvous_begin(r_earth, r_uranus); 
disp_var('thetaf', thetaf); 

WT = orbit.wait_time(thetaf, thetai, wt, wc, 1);
disp_var('WT', WT); 

WT = (thetaf - thetai)/(wt - wc); 
disp_var('WT', 365 - WT/86400); 


%% Problem 8 
ps = 1.225; %kg/m^3 
he = 122;
B = 0.1354; %km-1 

m = 3000;%kg
Ve = 7.7; %km/s
Cd = 1.7; 
S = 6.2/(1000*1000); % km^2 
l_over_d = 0.2; 
gamma_e = -1; % degrees 

% part b 
gamma = -6.7 * pi/180; 
 
V = sqrt(-2 * g_kms /(B*l_over_d * gamma));
disp_var('V', V); 

%% problem 10 

co2_day = 4*(8*0.62 + 2*7.18 + 6*7.18 + 8*1.04);
disp_var('co2', co2_day); 

lioh = co2_day * 0.75;

