%% Constants 
R = 6378; % km, Earth Radius
u = 3.986e5; % km3/s2
g = 9.81/1000; % km/s^2
orbit = orbital_mechanics_formulas();

function disp_var(var, value)
    disp([var, ' = ', num2str(value)]);
end

%% Problem 5
rp = 7500; % km;
T = 5 * 60 * 60; % hours to seconds
disp('Problem 5')

% Part a
a = (((T/(2*pi()))^2)*u)^(1/3);
disp_var('a', a)

% Part B 
e = 1 - rp/a;
disp_var('e', e);

% Part C
b = a * sqrt(1-e^2);
disp_var('b', b);

% Part D
p = a*(1-e^2);
disp_var('p', p);

% Part E 
delta_t = 3 * 60 * 60; % hours to seconds 
m = sqrt(u/(a^3))*delta_t; % mean anomaly
disp_var('m', m)

f = @(E) E - e * sin(E) - 3.77;
df = @(E) 1 - e * cos(E);
% Newton-Raphson iteration
E = 4;
max_iter = 100;
tol = 1e-3;
for i = 1:max_iter
    E_new = E - f(E) / df(E);  % Newton-Raphson formula
    if abs(E_new - E) < tol  % Convergence check
        break;
    end
    E = E_new;
end

disp_var('E', E);

f = 2*atan((sqrt(1+e)/sqrt(1-e)) * tan(E/2));
f_degrees = f * 180/pi;
disp_var('f', f);
disp_var('f_degrees', f_degrees);

disp('#############')


%% Problem 6
disp('Problem 6')
% Asteroid 1 
disp('Asteroid 1')
h = 1e5; % km^2/s
E = 1; % (km/s)^2

rp = roots([E, u, -(h^2)/2]);
disp('rp=')
disp(rp)

vinf = sqrt(2*E);
disp_var('vinf', vinf);
disp('#############')

% Asteroid 2
disp('Asteroid 2')
h = 1e5;
E = 100; 

rp = roots([E, u, -(h^2)/2]);
disp('rp = ')
disp(rp);

v_impact = sqrt(2*(E+u/R));
disp_var('v_impact', v_impact);
disp('#############')

% Asteroid 3 
disp('Asteroid 3')
E = 0;
h = 7e4;

rp = h^2/(2*u);
disp_var('rp', rp);

v_impact = sqrt(2*(E+u/R));
disp_var('v_impact', v_impact);
disp('#############')

% Asteroid 4 
disp('Asteroid 4')
E = 0;
h = 8e4;

rp = h^2/(2*u);
disp_var('rp', rp);

v_inf = sqrt(2*E);
disp_var('v_inf', v_inf);
disp('#############')

% Asteroid 5 
disp('Asteroid 5')
E = 10;
h = 8e4;

rp = roots([E, u, (-h^2/2)]);
disp('rp=');
disp(rp);

v_inf = sqrt(2*E);
disp_var('v_inf', v_inf);

disp('#############')

