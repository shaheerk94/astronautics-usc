R = 6378; % KM
u = 3.986e5; % km3/s2 

orbit = orbital_mechanics_formulas();


%% Problem 1
disp("##########################################")
disp("PROBLEM 1")
rp = R + (140 * 1.852); % rp = R + hp, 140 mi to km
disp(['rp = ', num2str(rp)])

ra = R + (150 * 1.852); % ra = R + ha, 150 mi to km
disp(['ra = ', num2str(ra)])

e = (ra-rp)/(ra+rp); 
disp(['e = ', num2str(e)])

p = rp*(1+e);
disp(['p = ', num2str(p)])

vp = orbit.perigee_velocity(e, u, p);
orbit.disp_var('vp', vp);

va1 = orbit.apogee_velocity(e, u ,p);
disp(['va = ', num2str(va)])
disp("##########################################")

%% Problem 2
disp("Problem 2")
P = 106 * 60; % 106mins to seconds 
rp = R + 200; % hp=200km
a = ((P*sqrt(u))/(2*pi))^(2/3); 
disp(['a = ', num2str(a)])

ra = 2*a -rp;
disp(['ra = ', num2str(ra)])

ha = ra - R; 
disp(['ha = ', num2str(ha)])
disp("##########################################")
%% Problem 3
disp("Problem 3")
a = 400 + R;
dv = 0.03; % m/s ->km/s
da = 2*sqrt((a^3)/u) * dv;
disp(['da = ', num2str(da)])
disp("##########################################")

%% Problem 4
disp("Problem 4")
r1 = 6570;
r2 = 42160; 
at = (r1 + r2)/2; 
dv = sqrt(u/r2) - sqrt(u*(2/r2-1/at)) + sqrt(u*(2/r1-1/at)) - sqrt(u/r1);
disp(['dv = ', num2str(dv)])

TOF = (pi*sqrt((at^3)/u))/3600; % hours
disp(['TOF = ', num2str(TOF)])
disp("##########################################")


%% Problem 5 
disp("Problem 5");
rp = R + (140 * 1.852); 
disp(['rp = ', num2str(rp)])
ra = R + (150 * 1.852);
disp(['ra = ', num2str(ra)])

ac = (ra + rp)/2;
disp(['a = ', num2str(a)])

at = R + (200 * 1.852);
disp(['at = ', num2str(at)])

dthetac = 2*pi*(1 - (ac/at)^(3/2));
disp(['dthetac = ', num2str(dthetac)])

dsc = dthetac * at;
disp(['dsc = ', num2str(dsc)])
disp("##########################################")

%% Problem 6 
ac = 6570;
at = 42160; 

wt = sqrt(u/(at^3));
disp(['wt = ', num2str(wt)])

wc = sqrt(u/(ac^3));
disp(['wc = ', num2str(wc)])

thetaf = pi * (1-sqrt(((ac + at)/(2*at))^3));
disp(['thetaf = ', num2str(thetaf)])

WT = (thetaf - pi)/(wt-wc);
disp(['WT = ', num2str(WT), ' seconds'])
disp(['WT = ', num2str(WT/60), ' minutes'])
