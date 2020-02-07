function [sol_2, FLAG_terminate] = ROB535_ControlsProject_part2_Team1 (TestTrack,Xobs_seen,curr_state)

    cline_interp = interp_cline(TestTrack.cline,2);
    
    centerLineTheta = TestTrack.theta;
    FLAG_terminate = 0;
    
    % Times in seconds
    total_time = .51;  % was .7
    dt = .01;
    
    gains = [5, 10000;
        0,0;
        0,0];
    % speed_ext = [10,25]; % speeds were 15 70
    %speed_ext = [10, 15]; % [11, 14]; % [10,10] worked %  [6,10] % [7,15]
    
    speed_ext = [11,14]; % 11,13
    lookahead = 2;
    centerline_strength = 0.1;

    % Gains
    kp1 = gains(1,1);
    ki1 = gains(2,1);
    kd1 = gains(3,1);
    kp2 = gains(1,2);
    ki2 = gains(2,2);
    kd2 = gains(3,2);
    
    
    % Sum Error
    sum_error = zeros(1,6);
    prev_error = zeros(1,6);
    U = [0, 0];

    n_elements = size(TestTrack.cline, 2);
%     % d_theta
    d_theta = 1./ abs(centerLineTheta(:, 2:end) - centerLineTheta(:, 1:end-1));
    d_theta(:,n_elements) = d_theta(:,n_elements - 1);
    d_theta_ext = [min(d_theta), max(d_theta)];
    
    
    % Number of time steps
    nsteps = total_time/dt;
    
    % Get intial line to track from obstacles
    %closest_obstacle_center = get_closest_obstacle(curr_state,Xobs_seen);
    closest_obstacle_center = get_next_obstacle(curr_state,Xobs_seen,cline_interp);
    centerLine = get_lane_from_obstacle_center( closest_obstacle_center, TestTrack);
    Cnormals = get_Cnormals(centerLine);
    
    for i=1:nsteps
        
        % Update previous values to avoid extra calculations
        last_closest_obs_center = closest_obstacle_center;
        last_centerLine = centerLine;
        
        % Get new closest obstacle and update lane and cnormals if
        % necessary
        %closest_obstacle_center = get_closest_obstacle(curr_state,Xobs_seen);
        closest_obstacle_center = get_next_obstacle(curr_state,Xobs_seen,cline_interp);
        if last_closest_obs_center ~= closest_obstacle_center
            centerLine = get_lane_from_obstacle_center( closest_obstacle_center, TestTrack);
            if last_centerLine ~= centerLine
                Cnormals = get_Cnormals(centerLine);
            end
        end

        
        nearest_goal_id = knnsearch(centerLine', [curr_state(1), curr_state(3)]);
        lookahead_goal_id = nearest_goal_id + lookahead;
        if lookahead_goal_id > size(centerLine, 2)
            lookahead_goal_id = size(centerLine, 2);
        end
        nearest_goal = centerLine(:, nearest_goal_id);
        length_scale = (d_theta(:, nearest_goal_id) - d_theta_ext(1))/(d_theta_ext(2) - d_theta_ext(1));
        goal_speed = length_scale * (speed_ext(2) - speed_ext(1)) + speed_ext(1);
        goal_state = [nearest_goal(1), goal_speed, nearest_goal(2), 0, 0, 0];
        error = goal_state - curr_state;
        goal_state(5) = centerLineTheta(:, lookahead_goal_id) -centerline_strength*dot(Cnormals(:, nearest_goal_id), [error(1);error(3)], 1); %atan2(error(3),error(1));
        error = goal_state - curr_state;
        sum_error = sum_error + error;
        d_error = (error - prev_error)/0.01;
        delta_f = kp1 * error(5) + ki1 * sum_error(5) + kd1 * d_error(5);
        Fx = kp2 * error(2) + ki2 * sum_error(2) + kd2 * d_error(2);
        if delta_f > 0.5
            delta_f = 0.5;
        end
        if delta_f < -0.5
            delta_f = -0.5;
        end
        %Solve for trajectory
        T=(i-1)*0.01:0.01:i*0.01;
        U = [U;[delta_f, Fx]];
        
        U_sim = U(end-1:end, :);
        if mod(i,1000) == 999
            [Y,~]=forwardIntegrateControlInput(U,x0);
        else
            [~,Y]=ode45(@(t,x)bike(t,x,T,U_sim),T,curr_state);
        end
        curr_state = Y(end,:);
        prev_error = error;
        if ((curr_state(1) > TestTrack.cline(1,end)) && (curr_state(3) > TestTrack.cline(2,end))) == 1
            FLAG_terminate = 1;
        end
    end
    
    sol_2 = U;

end



function dzdt=bike(t,x,T,U)
%constants
Nw=2;
f=0.01;
Iz=2667;
a=1.35;
b=1.45;
By=0.27;
Cy=1.2;
Dy=0.7;
Ey=-1.6;
Shy=0;
Svy=0;
m=1400;
g=9.806;

%generate input functions
%delta_f=interp1(T,U(:,1),t,'previous','extrap');
%F_x=interp1(T,U(:,2),t,'previous','extrap');
delta_f = U(end,1);
F_x = U(end,2);

%slip angle functions in degrees
a_f=rad2deg(delta_f-atan2(x(4)+a*x(6),x(2)));
a_r=rad2deg(-atan2((x(4)-b*x(6)),x(2)));

%Nonlinear Tire Dynamics
phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

F_zf=b/(a+b)*m*g;
F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

F_zr=a/(a+b)*m*g;
F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

F_total=sqrt((Nw*F_x)^2+(F_yr^2));
F_max=0.7*m*g;

if F_total>F_max
    
    F_x=F_max/F_total*F_x;
  
    F_yr=F_max/F_total*F_yr;
end

%vehicle dynamics
dzdt= [x(2)*cos(x(5))-x(4)*sin(x(5));...
          (-f*m*g+Nw*F_x-F_yf*sin(delta_f))/m+x(4)*x(6);...
          x(2)*sin(x(5))+x(4)*cos(x(5));...
          (F_yf*cos(delta_f)+F_yr)/m-x(2)*x(6);...
          x(6);...
          (F_yf*a*cos(delta_f)-F_yr*b)/Iz];
end


function Cnormals = get_Cnormals(centerLine)
    n_elements = size(centerLine, 2);
    projectors = (centerLine(:,2:n_elements) - centerLine(:,1:n_elements-1));
    projectors = projectors ./ vecnorm(projectors);
    Cnormals = [0, 1; -1, 0] * projectors;
    Cnormals(:,n_elements) = Cnormals(:,n_elements - 1);
end

function closest_obstacle_center = get_closest_obstacle(curr_state, Xobs_seen)
    num_cells = numel(Xobs_seen);
    if num_cells == 0
        closest_obstacle_center = [0;0];
        %disp("NO CELLS!");
    else
        min_dist = 999999;
        for i = 1:num_cells
            obstacle = Xobs_seen{i};
            obstacle_center = mean(obstacle)';
            obstacle_dist = norm(obstacle_center - [curr_state(1); curr_state(3)]);
            if obstacle_dist < min_dist
                min_dist = obstacle_dist;
                closest_obstacle_center = obstacle_center;
            end
        end
        
        %hold on;
        %scatter(closest_obstacle_center(1), closest_obstacle_center(2));
    end
end

function next_obstacle_center = get_next_obstacle(curr_state, Xobs_seen, cline_interp)
    num_cells = numel(Xobs_seen);
    next_obstacle_center = [0;0];
    if num_cells == 0
        %closest_obstacle_center = [0;0];
        %disp("NO CELLS!");
    else
        %centerLine = TestTrack.cline;
        centerLine = cline_interp;
        nearest_id_state = knnsearch(centerLine', [curr_state(1), curr_state(3)]);
        min_diff = 999999;
        for i = 1:num_cells
            obstacle = Xobs_seen{i};
            obstacle_center = mean(obstacle)';
            nearest_id_obs = knnsearch(centerLine', [obstacle_center(1), obstacle_center(2)]);
            id_diff = nearest_id_obs - nearest_id_state;
            if id_diff >= 0 && id_diff < min_diff
                next_obstacle_center = obstacle_center;
                min_diff = id_diff;
            end
        end
    end

end

function centerLine = get_lane_from_obstacle_center( closest_obstacle_center, TestTrack)
    if closest_obstacle_center == [0;0]
        %disp("Center")
        centerLine = TestTrack.cline;
    else
        
        n_elements = size(TestTrack.cline, 2);
        projectors = (TestTrack.cline(:,2:n_elements) - TestTrack.cline(:,1:n_elements-1));
        projectors = projectors ./ vecnorm(projectors);
        Cnormals = [0, 1; -1, 0] * projectors;
        Cnormals(:,n_elements) = Cnormals(:,n_elements - 1);
        
        % needs to take in centerLine and Cnormals
        % pos_obs needs to be row for knn to work
        cline_idx = knnsearch(TestTrack.cline', closest_obstacle_center');
        obs_vec = TestTrack.cline(:,cline_idx) - closest_obstacle_center;
        norm_vec = Cnormals(:,cline_idx);
        sign_check = dot(obs_vec, norm_vec);
        % positive = RHS
        % negative = LHS
        if sign_check >= 0
            %disp("Right")
            centerLine = (TestTrack.br + TestTrack.cline)./2;
        else
            %disp("Left")
            centerLine = (TestTrack.bl + TestTrack.cline)./2;
        end
    end
end

function cline_interp = interp_cline(cline, num_interpl)
    N = size(cline,2);
    step_size = 1/num_interpl;
    
    x_points = cline(1,:);
    y_points = cline(2,:);
    new_ind = 1:step_size:N;
    
    x_interp = interp1(x_points, new_ind);
    y_interp = interp1(y_points, new_ind);
    
    cline_interp = [x_interp;y_interp];
end
