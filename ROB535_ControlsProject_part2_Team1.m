function [sol_2, FLAG_terminate] = ROB535_ControlsProject_part2_Team1 (TestTrack,Xobs_seen,curr_state)
%     leftBound = TestTrack.bl;
%     rightBound = TestTrack.br;
%     centerLine = TestTrack.cline;
    
    
    %Xobs_seen
    
    closest_obstacle_center = get_closest_obstacle(curr_state,Xobs_seen);
    
    centerLine = get_lane_from_obstacle_center( closest_obstacle_center, TestTrack);
    
    
    
%     num_cells = numel(Xobs_seen);
%     if num_cells == 0
%         centerLine = TestTrack.cline;
%         disp("NO CELLS!");
%     else
%         min_dist = 999999;
%         %min_dist_obstacle_idx = 0;
%         closest_obstacle_center = [0;0];
%         for i = 1:num_cells
%             %bstacle_center = [0;0];
%             obstacle = Xobs_seen{i};
%             obstacle_center = mean(obstacle)
%             %for j = 1:4
%             %    obstacle_center = obstacle_center + obstacle(j,:);
%             %end
%             %obstacle_center = obstacle_center./4;
%             obstacle_dist = norm(obstacle_center - [curr_state(1); curr_state(3)]);
%             if obstacle_dist < min_dist
%                 min_dist = obstacle_dist;
%                 %min_dist_obstacle_idx = i;
%                 closest_obstacle_center = obstacle_center;
%             end
%         end
%         % Now we have closest obstacle index
%         closest_obstacle_center 
%         
%         hold on;
%         scatter(closest_obstacle_center(1), closest_obstacle_center(2))
%         
%     end
    
    
    
    %dist_obs_x = 
    %dist_obs_y = 
    %dot(Cnormals(:, nearest_obs_id), [dist_obs_x;dist_obs_y], 1);
    
    
    % To-do: Change centerline to avg between center line and proper
    % boundary. e.g. centerline = rightlane
    
    centerLineTheta = TestTrack.theta;
    FLAG_terminate = 0;
    
    % Times in seconds
    total_time = .7;
    dt = .01;
    
    gains = [5, 10000;
        0,0;
        0,0];
    speed_ext = [10,25]; % speeds were 15 70
    lookahead = 2;
    centerline_strength = 0.1;

    % Gains
    kp1 = gains(1,1);
    ki1 = gains(2,1);
    kd1 = gains(3,1);
    kp2 = gains(1,2);
    ki2 = gains(2,2);
    kd2 = gains(3,2);
    
    %speed = 6;
    
    % Get cnormals
    n_elements = size(centerLine, 2);
    projectors = (centerLine(:,2:n_elements) - centerLine(:,1:n_elements-1));
    projectors = projectors ./ vecnorm(projectors);
    Cnormals = [0, 1; -1, 0] * projectors;
    Cnormals(:,n_elements) = Cnormals(:,n_elements - 1);
    
    % Sum Error
    sum_error = zeros(1,6);
    prev_error = zeros(1,6);
    U = [0, 0];
    %Y_checking = [];
    simulation_divergence = zeros(1,6);
    
    % Cline Distances
    cline_distances = vecnorm(centerLine(:, 2:end) - centerLine(:, 1:end-1));
    cline_distances(:,n_elements) = cline_distances(:,n_elements - 1);
    d_ext = [min(cline_distances), max(cline_distances)];
    
    % d_theta
    d_theta = 1./ abs(centerLineTheta(:, 2:end) - centerLineTheta(:, 1:end-1));
    d_theta(:,n_elements) = d_theta(:,n_elements - 1);
    d_theta_ext = [min(d_theta), max(d_theta)];
    
    
    % Number of time steps
    nsteps = total_time/dt;
    
    % Initialize inputs
    %U = zeros(nsteps,2);
    
    for i=1:nsteps
        nearest_goal_id = knnsearch(centerLine', [curr_state(1), curr_state(3)]);
        lookahead_goal_id = nearest_goal_id + lookahead;
        if lookahead_goal_id > size(centerLine, 2)
            lookahead_goal_id = size(centerLine, 2);
        end
        %disp("Tracking point (out of 246) and iteration: ");
        %disp([nearest_goal_id, i]);
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
        %U(i,:) = [delta_f, Fx];
        
        U_sim = U(end-1:end, :);
        if mod(i,1000) == 999
            [Y,~]=forwardIntegrateControlInput(U,x0);
            [~,Y1]=ode45(@(t,x)bike(t,x,T,U_sim),T,curr_state);
            %simulation_divergence = Y(end,:)-Y1(end,:);
        else
            [~,Y]=ode45(@(t,x)bike(t,x,T,U_sim),T,curr_state);
        end
        curr_state = Y(end,:);
        %Y_checking = [Y_checking; [curr_state, goal_state, error, simulation_divergence]];
        prev_error = error;
        if close_to_end(curr_state, TestTrack)
            FLAG_terminate = 1;
            %last_input = U(end,:);
            last_input = [0, 4000];
            for i = 1:100
                U = [U;last_input];
            end
            sol_2 = U;
            break
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


function is_close = close_to_end(curr_state,TestTrack)
    x = curr_state(1);
    y = curr_state(3);
    last_point = TestTrack.cline(:,end);
    dist = norm(last_point-[x;y]);
    is_close = ( dist < 2 );
end

function closest_obstacle_center = get_closest_obstacle(curr_state, Xobs_seen)
    num_cells = numel(Xobs_seen);
    if num_cells == 0
        closest_obstacle_center = [0;0];
        %centerLine = TestTrack.cline;
        disp("NO CELLS!");
    else
        min_dist = 999999;
        %min_dist_obstacle_idx = 0;
        
        for i = 1:num_cells
            %bstacle_center = [0;0];
            obstacle = Xobs_seen{i};
            obstacle_center = mean(obstacle)';
            %for j = 1:4
            %    obstacle_center = obstacle_center + obstacle(j,:);
            %end
            %obstacle_center = obstacle_center./4;
            obstacle_dist = norm(obstacle_center - [curr_state(1); curr_state(3)]);
            if obstacle_dist < min_dist
                min_dist = obstacle_dist;
                %min_dist_obstacle_idx = i;
                closest_obstacle_center = obstacle_center;
            end
        end
        %hold on;
        %scatter(closest_obstacle_center(1), closest_obstacle_center(2));
    end
end

function centerLine = get_lane_from_obstacle_center( closest_obstacle_center, TestTrack)
    if closest_obstacle_center == [0;0]
        disp("Center")
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
        if sign_check >= 0
            disp("Right")
            centerLine = (TestTrack.br + TestTrack.cline)./2;
        else
            disp("Left")
            centerLine = (TestTrack.bl + TestTrack.cline)./2;
        end
        % positive = RHS
        % negative = LHS


        
        % for testing
        %centerLine = TestTrack.cline;
    end
end

% %function [sol_2, FLAG_terminate] = ROB535_ControlsProject_part2_Team1 (TestTrack,Xobs_seen,curr_state)
% 
% %%% where sol 1 is the vehicle?s trajectory.
% 
% %% sol 2 is a vector of control inputs that will be passed to forwardIntegrateControlInput. 
% %As mentioned later in 6.2, sol 2 must have enough control inputs with time step 0.01 sec- onds to forward integrate vehicle dynamics for the next 0.5 second.
% 
% Nw=2;
% f=0.01;
% Iz=2667;
% a=1.35;
% b=1.45;
% By=0.27;
% Cy=1.2;
% Dy=0.7;
% Ey=-1.6;
% Shy=0;
% Svy=0;
% m=1400;
% g=9.806;
% 
% 
% % (as an N � 2 vector where the first column is ? and the second column is Fx
% 
% %% linearize 
% 
% syms x1 x2 x3 x4 x5 x6 u1 u2
%  
%  
% a_f=rad2deg(u1-atan2(x4+a*x6,x2));
% a_r=rad2deg(-atan2((x4-b*x6),x2));
%  
%  
% F_zf=b/(a+b)*m*g;
% F_yf=F_zf*Dy*sin(Cy*atan(By*((1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy)))))+Svy;
%  
% F_zr=a/(a+b)*m*g;
% F_yr=F_zr*Dy*sin(Cy*atan(By*((1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy)))))+Svy;
%  
%  
% dzdx= [x2*cos(x5)-x4*sin(x5);...
%           (-f*m*g+Nw*u2-F_yf*sin(u1))/m+x4*x6;...
%           x2*sin(x5)+x4*cos(x5);...
%           (F_yf*cos(u1)+F_yr)/m-x2*x6;...
%           x6;...
%           (F_yf*a*cos(u1)-F_yr*b)/Iz];
%  
% A = jacobian(dzdx,[x1;x2;x3;x4;x5;x6]);
% B = jacobian(dzdx,[u1;u2]);
%  
% %% Decision Variables 
% % ith # of state : 6, input : 2 ;
% % horizon of 0.5 s & 0.01 s = 50 timestep
% Ndec = 6*51 + 2*5;
% npred = 50;
% 
% % Staates : z  &  Input : u 
% %z = zeros(6*51 + 2*5, 1 );
% %dzdt=bike(t,x,T,U)
% 
% %% Equality Constraints
% % quadprog
% % Euler integration  
% xsize = 6*51;
% zsize = Ndec;
% Aeq=zeros(xsize,zsize);
% Beq=zeros(xsize,1);
% 
% 
% for i = 1 : npred
%     Am = A(i);
%     Bm = B(i);
%     
%     Aeq(6*(i)+1 :6*(i)+6 , 6*(i-1)+1:6*(i-1)+6) = Am ;
%     Aeq(6*(i)+1:6*(i)+6, 6*(i-1)+1+6) = -ones(6,1);  
%     
%     Aeq(6*(i)+1:6*(i)+6 ,xsize+2*(i-1)+1:xsize+2*(i-1)+2) = Bm;
% end
% 
% Aeq(1:6,1:6) = eye(6);
% Beq(1:6) = reshape(curr_state,6,1);
% 
% 
% %%Boundary Constraints
% % 4-by-2 matrix 
%     Lb = -inf(Ndec,1);
%     Ub = inf(Ndec,1);
%     j=0;
%     for k = 1: npred
%         Lb(xsize+2*(k-1) + 1) = 0 - U(1,i+j);
%         Lb(xsize+2*(k-1) + 2) = -0.5 - U(2,i+j);
%  
%         Ub(xsize+2*(k-1) + 1) = 1 -U(1,i+j);
%         Ub(xsize+2*(k-1) + 2) = 0.5 -U(2,i+j);
%         j = j +1;
%     end
% 
% 
% %% quadatic penalty
% 
% 
% 
% 
% %% FLAG terminate is a binary flag set by the function to indicate when to stop the simulation. Setting this flag to 1 at a planning iteration indicates that, after integrating forward vehicle dynamics using the control inputs from sol 2 for 0.5 second, the vehicle would have reached the goal. If this flag is set to 0, then it implies that the car would not reach the goal after forward integrating for 0.5 second. Note that, in the event that FLAG terminate is never set to 1, the simulation would stop after 20 minutes so that we can process all the teams files in a timely manner.
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %
