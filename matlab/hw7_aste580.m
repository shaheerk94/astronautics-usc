function hw7_aste580

    %% Problem 1 
    r1_0 = [0, 0, 0];
    v1_0 = [0, 0, 0];
    r2_0 = [1, 0, 0];
    v2_0 = [0, sqrt(6), 0];
    r3_0 = [0.5, sqrt(3)/2, 0];
    v3_0 = [-3*sqrt(2)/2, sqrt(6)/2, 0];

    range = 0:pi/1000:pi;
    init_values = [r1_0 v1_0 r2_0 v2_0 r3_0 v3_0];
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'InitialStep', 1e-13);
    [t, drdt_dvdt] = ode45(@three_body_equations, range, init_values, options);

    r1 = [drdt_dvdt(:,1), drdt_dvdt(:,2), drdt_dvdt(:,3)];
    v1 = [drdt_dvdt(:,4), drdt_dvdt(:,5), drdt_dvdt(:,6)];

    r2 = [drdt_dvdt(:,7), drdt_dvdt(:,8), drdt_dvdt(:,9)];
    v2 = [drdt_dvdt(:,10), drdt_dvdt(:,11), drdt_dvdt(:,12)];
    
    r3 = [drdt_dvdt(:,13), drdt_dvdt(:,14), drdt_dvdt(:,15)];
    v3 = [drdt_dvdt(:,16), drdt_dvdt(:,17), drdt_dvdt(:,18)];

    r_cm = (r1 + 2.*r2 + 3.*r3)./6;
    % Indices where t = 0, pi/4, pi/2, 3pi/4
    t_vals = [0, pi/4, pi/2, 3*pi/4];
    t_indices = find(ismember(t, t_vals));
    
    % Relative positions: r_i - r_cm for r1, r2, r3
    r_all = {r1, r2, r3};
    r_rel = cell(3,1);
    for j = 1:3
        r_rel{j} = r_all{j}(t_indices,:) - r_cm(t_indices,:);
    end
    
    % Construct closed-loop polygons (each row = one time)
    polys_x = cell(4,1);
    polys_y = cell(4,1);
    for k = 1:4
        px = zeros(1,4);
        py = zeros(1,4);
        for j = 1:3
            px(j) = r_rel{j}(k,1);
            py(j) = r_rel{j}(k,2);
        end
        px(4) = px(1);
        py(4) = py(1);
        polys_x{k} = px;
        polys_y{k} = py;
    end
    
    % Plot full relative trajectories
    figure(1); hold on;
    for j = 1:3
        plot(r_cm(:,1) - r_all{j}(:,1), r_cm(:,2) - r_all{j}(:,2));
    end
    
    % Plot closed-loop triangles at specific time steps
    for k = 1:4
        plot(polys_x{k}, polys_y{k}, 'LineWidth', 1);
    end
    grid on; 
    legend('m1', 'm2', 'm3', 't=0', 't=pi/4', 't=pi/2', 't=3pi/4')
    xlabel('x position')
    ylabel('y position')

    %% Problem 2 
    u = 1e-5:(0.5-1e-5)/10:0.5;
    % Define the potential derivative at y = 0
    dVdx = @(x, mu) x - (1 - mu) * (x + mu) ./ ((x + mu).^2).^(3/2) - mu * (x - 1 + mu) ./ ((x - 1 + mu).^2).^(3/2);
    
    L1 = zeros(size(u));
    L2 = zeros(size(u));
    L3 = zeros(size(u));
    
    for i = 1:length(u)
        mu = u(i);
        f = @(x) dVdx(x, mu);
        
        % Find L1, L2, L3 near their expected locations
        L1(i) = fzero(f, 1 - mu/3);     % Between two bodies
        L2(i) = fzero(f, 1.2);          % Beyond secondary
        L3(i) = fzero(f, -1);           % Beyond primary
    end
    
    % Plotting
    figure;
    plot(u, L1, 'r-', 'DisplayName', 'L1', 'LineWidth', 1.5);
    hold on;
    plot(u, L2, 'b-', 'DisplayName', 'L2', 'LineWidth', 1.5);
    plot(u, L3, 'k-', 'DisplayName', 'L3', 'LineWidth', 2.5);
    xlabel('\mu');
    ylabel('x-coordinate');
    legend;
    title('x-location of Lagrange points vs. \mu');
    grid on;


    
end

function drdt_dvdt = three_body_equations(~, y)
    r1 = y(1:3);
    v1 = y(4:6);
    r2 = y(7:9);
    v2 = y(10:12);
    r3 = y(13:15);
    v3 = y(16:18);

    r21 = r1 - r2;
    r31 = r1 - r3;
    r12 = r2 - r1;
    r32 = r2 - r3;
    r13 = r3 - r1;
    r23 = r3 - r2;

    % r1dot = v1
    drdt_dvdt(1, 1) = v1(1);
    drdt_dvdt(2, 1) = v1(2);
    drdt_dvdt(3, 1) = v1(3);

    % v1dot
    drdt_dvdt(4, 1) = -2/(norm(r21)^3) * r21(1) - 3/(norm(r31)^3) * r31(1);
    drdt_dvdt(5, 1) = -2/(norm(r21)^3) * r21(2) - 3/(norm(r31)^3) * r31(2);
    drdt_dvdt(6, 1) = -2/(norm(r21)^3) * r21(3) - 3/(norm(r31)^3) * r31(3);
    
    % r2dot = v2 
    drdt_dvdt(7, 1) = v2(1);
    drdt_dvdt(8, 1) = v2(2);
    drdt_dvdt(9, 1) = v2(3);
    
    % v2dot 
    drdt_dvdt(10, 1) = -1/(norm(r12)^3) * r12(1) - 3/(norm(r32)^3) * r32(1);
    drdt_dvdt(11, 1) = -1/(norm(r12)^3) * r12(2) - 3/(norm(r32)^3) * r32(2);
    drdt_dvdt(12, 1) = -1/(norm(r12)^3) * r12(3) - 3/(norm(r32)^3) * r32(3);

    % r3dot = v3
    drdt_dvdt(13, 1) = v3(1);
    drdt_dvdt(14, 1) = v3(2);
    drdt_dvdt(15, 1) = v3(3);

    % v3dot 
    drdt_dvdt(16, 1) = -1/(norm(r13)^3) * r13(1) - 2/(norm(r23)^3) * r23(1);
    drdt_dvdt(17, 1) = -1/(norm(r13)^3) * r13(2) - 2/(norm(r23)^3) * r23(2);
    drdt_dvdt(18, 1) = -1/(norm(r13)^3) * r13(3) - 2/(norm(r23)^3) * r23(3);
end