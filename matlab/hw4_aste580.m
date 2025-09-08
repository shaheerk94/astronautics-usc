function hw4_aste580

    % Problem 1 
    % Constants
    u = 398600; % km^3/s^2
    
    % Initial conditions
    r0 = [10000; 1000; 0]; % km
    v0 = [1; (sqrt(u / 10000) + 1); 1]; % km/s

    E = norm(v0)^2/2 - u/norm(r0);
    a = -u/(2*E);
    disp_var('a', a);
    
    % Time settings
    T = 2 * pi * sqrt(a^3 / u); % Orbital period
    tspan = [0, 10 * T];
    dt = 10 * 60; % 10 minutes in seconds
    
    % ODE solver
    y0 = [r0; v0];
    options = odeset('RelTol', 1e-11, 'AbsTol', 1e-12, 'Refine', 10);
    [t, y] = ode45(@(t, y) two_body_equations(t, y, u), 0:dt:tspan(2), y0, options);
    
    % Extract results
    r = y(:, 1:3);
    v = y(:, 4:6);
    
    % Energy and Angular Momentum
    energy = 0.5 * sum(v.^2, 2) - u ./ vecnorm(r, 2, 2);
    angular_momentum = cross(r, v, 2);
    
    % Plotting
    figure;
    subplot(2,2,1);
    plot(r(:,1), r(:,2)); grid on;
    xlabel('x (km)'); ylabel('y (km)'); title('Orbit in x-y plane');
    
    subplot(2,2,2);
    plot(r(:,1), r(:,3)); grid on;
    xlabel('x (km)'); ylabel('z (km)'); title('Orbit in x-z plane');
    
    subplot(2,2,3);
    plot(t/3600, energy); grid on;
    xlabel('Time (hours)'); ylabel('Energy (km^2/s^2)'); title('Energy vs Time');
    
    subplot(2,2,4);
    plot(t/3600, vecnorm(angular_momentum, 2, 2)); grid on;
    xlabel('Time (hours)'); ylabel('Angular Momentum (km^2/s)'); title('Angular Momentum vs Time');

    saveas(gcf, 'hw4_aste580_.png')

    figure;
    plot3(r(:,1), r(:,2), r(:,3), 'LineWidth', 1.5);
    grid on; xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
    title('3D Orbit Trajectory');

    %% Problem 2
    % Part D 
    u = 4e5; % km^3/s^2
    r1 = 600 + 6400; % km
    
    p = [2, 5, inf];
    di = [60, 45, 30];

    % Assuming r1 and u are defined earlier in your code
    for i = 1:length(di)
        for j = 1:length(p)
            % Calculate delta_v and delta_vt for each combination of p and di
            deltav = one_impulse_cost(r1, u, di(i));
            deltavt = three_impulse_cost(r1, u, p(j), di(i));

            if (p(j) == inf)
                deltavt = three_impulse_cost_inf_p(u, r1);
            end
        
            % Calculate the ratio R
            R = deltavt / deltav;
        
            % Print the results nicely
            fprintf('For p = %g, di = %d:\n', p(j), di(i));
            fprintf('  delta_v = %.4f, delta_vt = %.4f, R = %.4f\n\n', deltav, deltavt, R);
        end
    end

    %% Part E 
    i = 2 * acsc(1/(sqrt(2) -1));
    disp_var('i', i*180/pi);

    %% Part F 
    i = 2 * acsc(3) * 180/pi;
    disp_var('i', i);
    

end

function dydt = two_body_equations(~, y, u)
    r = y(1:3);
    v = y(4:6);
    r_mag = norm(r);
    a = -u / r_mag^3 * r;
    dydt = [v; a];
end

function deltav = one_impulse_cost(r, u, i)
    deltav = 2 * sqrt(u/r) * sind(i/2);
end

function deltaVT = three_impulse_cost(r, u, p, i) 
    deltaVT = 2 * sqrt(u/r) * ( sqrt(2*p/(1+p)) -1 + sind(i/2) * sqrt(2/(p*(1+p))) );
end

function deltaVT = three_impulse_cost_inf_p(u, r)
    deltaVT = 2 * sqrt(u/r) * (sqrt(2) -1);
end


function disp_var(var, value)
    disp([var, ' = ', num2str(value)]);
end

