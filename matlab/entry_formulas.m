%% General Entry Functions  
    % Definitions: 
    % a = acceleration 
    % Cl=coefficient of lift 
    % Cd=coefficient of drag
    % gamma = flight path angle
    % gamme_e = flight path angle @ entry interface 
    % h=altitude 
    % m=spacecraft mass 
    % p=atmospheric density
    % ps = surface atmospheric density 
    % S=lift surface area
    % Vr=velocity relative to atmosphere
    % Ve = velocity @ entry interface  
classdef entry_formulas
    methods(Static)
        function disp_var(var, value)
            disp([var, ' = ', num2str(value)]);
        end
        
        
        function L = lift_from_Cl(p, Vr, S, Cl)
            L = 0.5 * p * Vr^2 * S * Cl;
        end
        
        function D = drag_from_Cd(p, Vr, S, Cd)
            D = 0.5 * p * Vr^2 * S * Cd;
        end
        
        function ph = atmospheric_density(ps, B, h)
            ph = ps * log(-B * h);
        end
        
        
        %% Ballistic Entry 
        % instantaneous V as a function of inst. altitude, density, and
        % constants 
        % note gamma_e in degrees 
        function V = inst_ballistic_velocity(Ve, B, ps, gamma_e, S, Cd, m, h)
            V = Ve * log((-1 * ps * S * Cd * log(-B * h))/(2*B * sind(gamma_e) * m));
        end
        
        
        function amax = max_accel_ballistic(B, Ve, gamma_e, e)
        % max acceleration for ballistic entry 
        % note gamma_e in degrees 
            amax = (-1 * B * Ve^2)/(2 * e) * sind(gamma_e);
        end
        
        % altitude at which amax occurs 
        % note gamma_e in degrees 
        function hcrit = altitude_of_max_accel(B, S, Cd, m, ps, gamma_e)
            hcrit = (1/B) * log((-1 * S * Cd * ps)/(B * m * sind(gamma_e)));
        end
        
        % velocity at hcrit
        % note gamme_e in degrees 
        function vcrit = velocity_at_hcrit(Ve)
            vcrit = Ve/sqrt(e);
        end
        
        
        %% Glided Entry 
        
        % flight path angle as a function of inst. velocity
        function gamma = inst_glided_flight_path_angle(B, l_over_d, V, g)
            gamma = -2 / (B * l_over_d * (V^2/g));
        end
        
        % inst. velocity as a function of density and D=drag
        function V = inst_glided_velocity(g, R, S, Cd, m, l_over_d, p)
            V = sqrt((g*R)/(1 + (R*S*Cd*l_over_d*p)/(2*m)));
        end
        
        % inst. acceleration as a function of inst. V and constants 
        function a_g = inst_glided_accel(V, g, R, l_over_d)
            a_g = ((V^2/(g*R)) - 1)/l_over_d;
        end
        
        % total range as a function of inst velocity at entry interface and
        % constants 
        function s = total_glided_range(R, l_over_d, Ve, g)
            s = 0.5 * R * l_over_d * log(1/(1-(Ve^2)/(g*R)));
        end
        
        % cross range, for gliding vehicle only
        % note theta in degrees
        function s_cross = cross_ranging_dist_glided(R, l_over_d, theta)
            s_cross = R * pi^2 * (l_over_d^2)/48 * sind(2*theta);
        end
    end
end
