%% Problem 1 
format shortE
J2 = 0.452e-5;
uvenus = 324820; 
rvenus = 6052.75; % km
T = 225 * 3600 * 24; % seconds 

nv = 2*pi/T; 
disp(nv)

a = (1.5 * sqrt(uvenus) * rvenus^2 * J2/nv)^(2/7)
disp(a) 