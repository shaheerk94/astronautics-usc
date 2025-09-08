%% Constants 
R = 6378; % km, Earth Radius
u = 3.986e5; % km3/s2
g = 9.81/1000; % km/s^2
orbit = orbital_mechanics_formulas();

function disp_var(var, value)
    disp([var, ' = ', num2str(value)]);
end

%% Problem 4
rp = 640 + R;
vp = 10.4; % km/s
e = vp^2*rp/u -1;
disp_var('e', e);
a = rp/(1-e);
disp_var('a', a);
ra = a * (1+e);
disp_var('ra', ra);

p = 13361;
e = 0.9044;
a = p/(1-e^2);
disp_var('a', a)
T = 2*pi*sqrt(a^3/u);
disp_var('T', T);
T_hours = T/3600; 
disp_var('T_hours', T_hours);