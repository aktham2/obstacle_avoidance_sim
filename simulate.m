function [data, params] = simulate
    
    % Parameters
    % - gravity
    params.g = 9.81;
    % - mass
    params.m = 0.7;
    % - moment of inertia
    params.J = diag([0.004, 0.005, 0.007]);
    % - spar length
    params.l = 0.17;
    % - aerodynamic force and moment coefficients
    params.kF = 1e-5;
    params.kM = 1e-7;
    % - maximum spin rate
    params.sigmamax = (1e3);
    % - radius of bounding volume (sphere) around quadrotor
    params.r = 3*params.l;
    % - sample time
    params.dt = (1/50);
    
    % Simulation
    % - initial time
    t0 = 0;
    % - initial state
    o0 = [-1.5; 1; -2];
    theta0 = [0; 0; 0];
    v0 = [0; 0; 0];
    w0 = [0; 0; 0];
    x0 = [o0; theta0; v0; w0];
    % - final time
    t1 = 10;
    
    % Problem
    % - desired position
    o_desired = [-1.5; 1.5; -2];
    % - goal position
    o_goal = [1.5; -1.0; -1.5];
    % - obstacles
    %   * create an empty cell array to hold obstacles
    obst = {};
    %   * uncomment this line to add a single spherical obstacle (center, radius)
%      obst = AddObstacle_Sphere(obst, [0; 1; -2], 0.5)
    %   add n randomly located spherical obstacles:
    n = 20;
    obst = AddObstacle_RandomSpheres(obst, n, 0.1, 0.5, 2.5, 2.5, 2.5, ...
                                     o_desired, o_goal, params);
	
    % DESIGN
    % - controller    
    params.xe = [o_desired;zeros(9,1)];
    params.ue = [0;0;0;params.m*params.g];
    [Ad,Bd] = getAdBd(params.xe,params.ue,params.dt,params.m,params.g,params.J);
    Q = diag([100,100,100,1,1,1,1,1,1,1,1,1]);
    R = diag([1,1,1,0.25]);
    params.K = dlqr(Ad,Bd,Q,R);
    % - planner
    params.k_att = 1;
    params.b_att = 1;
    params.k_rep = 1.1;
    params.b_rep = .1;
    params.k_des = 1*params.dt;
    params.b_des = 1*params.dt;
    
    % Data to store
    data.t = [t0];
    data.x = [x0];
    data.u = [];
    data.o_desired = [o_desired];
    data.o_goal = [o_goal];
    data.obst = {obst};
    
    % Iterate over t1/dt sample intervals.
    for i = 1:round(t1/params.dt)

        % Current time and state
        t = data.t(:, end);
        x = data.x(:, end);
        
        % Planner (update o_desired)
        o_desired = planner(t, x, o_desired, o_goal, obst, params);
        params.xe = [o_desired; 0; 0; 0; 0; 0; 0; 0; 0; 0]; 
        
        % Controller (update u)
        u = controller(t, x, o_desired, params);
        
        % Simulator (update t and x)
        [tsol, xsol] = ode45(@(t, x) h(t, x, u, params.g, params.m, params.J), [t t+params.dt], x);
        
        % Store data
        data.o_desired(:, end+1) = o_desired;
        data.u(:, end+1) = u;
        data.t(:, end+1) = tsol(end, :)';
        data.x(:, end+1) = xsol(end, :)';
        data.obst{:, end+1} = obst;
    end

    % Get position and orientation
    data.o = data.x(1:3, :);
    data.theta = data.x(4:6, :);

end

function obst = AddObstacle_Sphere(obst, p, s)
    obst{end+1} = struct('type', 1, 'p', p, 's', s);
end

function obst = ...
    AddObstacle_RandomSpheres(obst, ...
                              n, smin, smax, dx, dy, dz, ...
                              o_start, o_goal, params)

	% n = number of spheres to add
    % smin = minimum radius
    % smax = maximum radius
    % o_start = start position
    % o_goal = goal position
    % (-dx, dx), (-dy, dy), (-dz, 0) = room dimensions

    % Loop to add n spheres one at a time
    for i=1:n

        % Keep trying until we add a sphere that does not interfere
        % with the start or goal position of the quadrotor.
        while(1)
            % Random center position in the room
            p = [dx;dy;dz].*([-1;-1;-1]+[2;2;1].*rand(3,1));
            % Random radius
            s = smin+(smax-smin)*rand;
            % Check for non-interference
            ds = norm(p-o_start)-(params.r+s);
            dg = norm(p-o_goal)-(params.r+s);
            if ((ds>0)&&(dg>0.5))
                obst = AddObstacle_Sphere(obst, p, s);
                break;
            end
        end

    end
end

function xdot = h(t, x, u, g, m, J)
    
    %%%%%%%%%
    % FIXME %
    %%%%%%%%%
       c1 = cos(x(4));
    s1 = sin(x(4));
    c2 = cos(x(5));
    s2 = sin(x(5));
    c3 = cos(x(6));
    s3 = sin(x(6));
    odot = [x(7);x(8);x(9)];
    N = [-s2 0 1;
        c2*s3 c3 0;
        c2*c3 -s3 0];
    tdot = N\[x(10);x(11);x(12)];
    Rz = [c1 -s1 0; s1 c1 0; 0 0 1];
    Ry = [c2 0 s2; 0 1 0; -s2 0 c2];
    Rx = [1 0 0; 0 c3 -s3; 0 s3 c3];
    R = Rz*Ry*Rx;
    vdot = [0;0;g] - R*[0;0;u(4)/m];
    wdot = J\([u(1); u(2); u(3)]-[0 -x(12) x(11); x(12) 0 -x(10); -x(11) x(10) 0]*J*[x(10);x(11);x(12)]);
    xdot = [odot; tdot; vdot; wdot];
    
end

function u = controller(t, x, o_desired, params)
    
    %%%%%%%%%
    % FIXME %
    %%%%%%%%%
    
    u = -params.K*(x-params.xe) + params.ue;    
    W = [params.l*params.kF -params.l*params.kF 0 0;
        0   0  params.l*params.kF -params.l*params.kF;
        params.kM params.kM -params.kM -params.kM;
        params.kF params.kF params.kF params.kF];
    s = W\u;
    
     for i=1:length(s)
        if s(i)>params.sigmamax^2
            s(i) = params.sigmamax^2;
        end
        if s(i)<0
            s(i)=0;
        end
     end
     u = W*s;
end

function o_desired = planner(t, x, o_desired, o_goal, obst, params)
    
    %%%%%%%%%
    % FIXME %
    %%%%%%%%%
    
    q = o_desired;
    q_goal = o_goal;
    r = params.r;
    
    
    if(norm(q-q_goal) <= params.b_att)
        gradf = params.k_att * (q-q_goal);
    else
        gradf = params.b_att*params.k_att* (q-q_goal)/norm(q-q_goal);
    end
    for i=1:length(obst)
            d = norm(q-obst{i}.p) - (r+obst{i}.s);
        if( d <= params.b_rep)
            gradf = gradf + -params.k_rep * (1/d - 1/params.b_rep) * (1/d)^2 ...
                     * ((q-obst{i}.p)/norm(q-obst{i}.p));
        end
    end
	
	
    if(norm((params.k_des*gradf)) <= params.b_des)
        q = q - params.k_des * gradf;
    else
        q = q - params.b_des*(gradf/norm(gradf));
    end
    o_desired = q;
    % 

end







