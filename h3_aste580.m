%% Constants 
R = 6378; % km, Earth Radius
u = 3.986e5; % km3/s2
g = 9.81/1000; % km/s^2
orbit = orbital_mechanics_formulas();

function disp_var(var, value)
    disp([var, ' = ', num2str(value)]);
end

%% Problem 1 
f = acos(2/32) -1;
f_deg = f * 180/pi;
disp_var('f_deg', f_deg);
v = sqrt(1/2 * (1 + 2*cos(f) +1));
disp_var('v', v); 

%% Problem 2 
r_t0 = [10000, 1000, 0]; % km
v_t0 = [1, sqrt(u/10000) + 1, 1]; % km/s

[a, e, i, Energy, h, Omega, w, f, E, M, n, Period] = orbital_elements(r_t0, v_t0, u);
disp_var('Semi-Major Axis', a); 
disp_var('Eccentricity', e); 
disp_var('Inclination', i);
disp_var('RAAN', Omega);
disp_var('Argument of Periapsis', w);
disp_var('Mean Anomaly', M);
disp_var('Orbital Energy', Energy);
disp_var('Angular Momentum Vector', h);
disp_var('Angular Momentum Magnitude', h_mag);
disp_var('Mean Motion', n);
disp_var('Period', Period);
disp_var('True Anomaly', f);
disp_var('Eccentric Anomaly', E);

disp_var('Period', Period);

function [a, e, i, Energy, h, Omega, w, f, E, M, n, Period] = orbital_elements(r_t0, v_t0, u)
    r_t0_mag = norm(r_t0); 
    v_t0_mag = norm(v_t0); 
    Energy = v_t0_mag^2/2 - u/r_t0_mag;
    h = cross(r_t0, v_t0); 
    h_mag = norm(h);
    p = h_mag^2/u;
    a = -1* u/(2*Energy);
    e = sqrt(1-p/a); 
    e_vector = (1/u) * cross(v_t0, h) - r_t0/r_t0_mag;
    i = acos(h(3)/h_mag);
    N = cross([0, 0, 1], h);
    n_mag = norm(N); 
    if n_mag ~= 0
        if N(2) >= 0
            Omega = acos(N(1)/n_mag);
        else
            Omega = 2 * pi - acos(N(1)/n_mag);
        end
    else
        Omega = 0;
    end
    if e_vector(3) >= 0
        w = acos(dot(N, e_vector)/m(n_mag * e));
    else
        w = 2 * pi - acos(dot(N, e_vector)/(n_mag * e));
    end

    f = acos((p/r_t0_mag -1)/e);
    E = 2 * atan((sqrt(1-e)/sqrt(1+e)) * tan(f/2));
    M = E - e*sin(E);
    n = sqrt(u/a^3);
    Period = 2 * pi * sqrt(a^3/u);
end