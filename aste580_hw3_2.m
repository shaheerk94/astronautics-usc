% MATLAB script for computing orbital elements

% Given initial conditions
r0 = [10000, 1000, 0]; % Position vector in km
v0 = [1, sqrt(398600.4418 / 1000) + 1, 1]; % Velocity vector in km/s
mu = 398600.4418; 

% Compute specific angular momentum
h = cross(r0, v0);
h_mag = norm(h);

% Compute eccentricity vector
e_vec = (cross(v0, h) / mu) - (r0 / norm(r0));
e = norm(e_vec);

% Compute semi-major axis
a = 1 / ((2 / norm(r0)) - (norm(v0) ^ 2 / mu));

% Compute inclination
i = acos(h(3) / h_mag);

% Compute right ascension of ascending node (RAAN)
n = cross([0, 0, 1], h);
n_mag = norm(n);
if n_mag ~= 0
    Omega = acos(n(1) / n_mag);
    if n(2) < 0
        Omega = 2 * pi - Omega;
    end
else
    Omega = 0;
end

% Compute argument of periapsis
if e ~= 0 && n_mag ~= 0
    omega = acos(dot(n, e_vec) / (n_mag * e));
    if e_vec(3) < 0
        omega = 2 * pi - omega;
    end
else
    omega = 0;
end

% Compute true anomaly
if e ~= 0 && dot(e_vec, r0) / (e * norm(r0)) <= 1
    nu = acos(dot(e_vec, r0) / (e * norm(r0)));
    if dot(r0, v0) < 0
        nu = 2 * pi - nu;
    end
else
    nu = 0;
end

% Compute eccentric anomaly
if e ~= 0 && (1 + e) / (1 - e) > 0
    E = 2 * atan(tan(nu / 2) / sqrt((1 + e) / (1 - e)));
else
    E = nu;
end

% Compute mean anomaly
M = E - e * sin(E);

% Compute energy
energy = (norm(v0) ^ 2 / 2) - (mu / norm(r0));

% Compute mean motion
n = sqrt(mu / a ^ 3);

% Compute period
T = 2 * pi / n;

% Display results
fprintf('Semi-major axis (a): %.6f km\n', a);
fprintf('Eccentricity (e): %.6f\n', e);
fprintf('Inclination (i): %.6f degrees\n', rad2deg(i));
fprintf('RAAN (Ω): %.6f degrees\n', rad2deg(Omega));
fprintf('Argument of periapsis (ω): %.6f degrees\n', rad2deg(omega));
fprintf('Mean anomaly (M): %.6f degrees\n', rad2deg(M));
fprintf('Specific orbital energy: %.6f km^2/s^2\n', energy);
fprintf('Specific angular momentum: %.6f km^2/s\n', h_mag);
fprintf('Mean motion (n): %.6f rad/s\n', n);
fprintf('Orbital period (T): %.6f seconds\n', T);
fprintf('True anomaly (ν): %.6f degrees\n', rad2deg(nu));
fprintf('Eccentric anomaly (E): %.6f degrees\n', rad2deg(E));
