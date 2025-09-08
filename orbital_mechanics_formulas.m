classdef orbital_mechanics_formulas 
    methods(Static) 

        function disp_var(var, value)
            disp([var, ' = ', num2str(value)]);
        end

        function km = convert_naut_miles(miles)
            km = miles * 1.852;
        end

        function km = convert_miles(miles)
            km = miles * 1.60934;
        end

        % change in semi-major axis from a circular orbit given delta v
        function delta_a = circular_delta_a(a, u, delta_v)
            delta_a = 2 * sqrt(a^3/u) * delta_v;
        end

        % orbital period given semimajor axis 
        function P = orbital_period(u, a)
            P = 2 * pi/sqrt(u) * a^(3/2);
        end
        
        % orbital velocity as a function of inst. r and a 
        function v = orbital_velocity(u, r, a)
            v = sqrt(u * (2/r - 1/a));
        end

        % perigee velocity 
        
        function vp = perigee_velocity(e, u, p)
            vp = (1 + e) * (sqrt(u/p));
        end

        % apogee velocity 
        function va = apogee_velocity(e, u, p)
            va = (1-e) * sqrt(u/p);
        end
        
        % first hohmann transfer dv, from circular 1 to elliptical 
        function dv1_hohmann = first_hohmann_dv(u, r1, a)
            dv1_hohmann = sqrt(u * (2/r1 - 1/a));
        end

        % second hohmann transfer dv, from elliptical to circular 2 
        function dv2_hohmann = second_hohmann_dv(u, r2, a)
            dv2_hohmann = sqrt(u * (2/r2 - 1/a)) - sqrt(u/r2);
        end

        % total hohmann transfer dv
        function dv_tot_hohman = total_hohman_dv(u, r1, r2, at)
            dv_tot_hohman = sqrt(u/r2) - sqrt(u*(2/r2 - 1/at)) + ...
            sqrt(u*(2/r1 - 1/at)) - sqrt(u/r1);
        end

        function tof = time_of_flight_hohmann(at, u)
            tof = pi * sqrt((at^3)/u);
        end

        % dv for inclination change 
        function dv_inclination = dv_inclination_change(v, thetad)
            dv_inclination = 2 * v * sind(thetad/2);
        end

        % change in angle between chaser and target per rev of chaser
        function dtheta_chaser = dtheta_chaser(ac, at)
            dtheta_chaser = 2 * pi * (1-((ac/at)^(3/2)));
        end

        % change in arc length across target orbit per rev of chaser 
        function ds_chaser = arc_length_from_phase_change(ac, at)
            ds_chaser = 2 * pi * at * (1 - ((ac/at)^(3/2)));
        end
        
        % angular speed given a (plug in ac for chaser, at for target) 
        function w = angular_speed(u, a)
            w = sqrt(u/(a^3));
        end

        % angle @ which to begin transfer for rendezvous
        function thetaf = phase_angle_rendezvous_begin(ac, at)
            thetaf = pi * (1 - ((ac + at)/(2*at))^(3/2));
        end

        % wait time given initial angles and speeds 
        function WT = wait_time(thetaf, thetai, wt, wc, k)
            WT = (thetaf - thetai)/(wt - wc) + (2 * pi * k)/(wt -wc);
        end

        % total time of flight in transfer orbit 
        function TOF = time_of_flight_rendezvous(ac, at, u)
            atrans = (ac + at)/2;
            TOF = pi * sqrt(atrans^3/u);
        end

        % semimajor axis of chaser orbit in co-planar rendezvous
        function ac = ac_coplanar(u, thetai, wt) 
            ac = (u * ((thetai/(2*pi*wt))^2))^(2/3);
        end

        % first dv required for coplanar rendezvous
        function dv = dv_initial_coplanaar(u, at, ac)
            dv = sqrt(u * (2/at - 1/ac)) - sqrt(u/at);
        end

        % total DV is just the double of the initial DV 
        function dv_tot = dv_total_coplanar(u, ac, at)
            dv_tot = 2 * abs(sqrt(u * (2/at - 1/ac)) - sqrt(u/at));
        end

    end
end