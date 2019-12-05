function [Y_checking, U] = PDControlPart1(nsteps, x0, gains, speed_ext, lookahead, centerline_strength)
    load TestTrack.mat
    % Track boundaries
    leftBound = TestTrack.bl;
    rightBound = TestTrack.br;
    centerLine = TestTrack.cline;
    centerLineTheta = TestTrack.theta;
    %{
    for i=1:subdivision_num
        centerLine = subdivideTrack(centerLine);
        centerLineTheta = subdivideTrack(centerLineTheta);
    end
    %}
    
    kp1 = gains(1,1);
    ki1 = gains(2,1);
    kd1 = gains(3,1);
    kp2 = gains(1,2);
    ki2 = gains(2,2);
    kd2 = gains(3,2);
    
    % Required timestep
    setGlobaldt(0.01);

    curr_state = x0;
    n_elements = size(centerLine, 2);
    projectors = (centerLine(:,2:n_elements) - centerLine(:,1:n_elements-1));
    projectors = projectors ./ vecnorm(projectors);
    Cnormals = [0, 1; -1, 0] * projectors;
    Cnormals(:,n_elements) = Cnormals(:,n_elements - 1);

    sum_error = zeros(1,6);
    prev_error = zeros(1,6);
    U = [0, 0];
    Y_checking = [];
    simulation_divergence = zeros(1,6);
    
    cline_distances = vecnorm(centerLine(:, 2:end) - centerLine(:, 1:end-1));
    cline_distances(:,n_elements) = cline_distances(:,n_elements - 1);
    d_ext = [min(cline_distances), max(cline_distances)];
    
    d_theta = 1./ abs(centerLineTheta(:, 2:end) - centerLineTheta(:, 1:end-1));
    d_theta(:,n_elements) = d_theta(:,n_elements - 1);
    d_theta_ext = [min(d_theta), max(d_theta)];

    for i=1:nsteps
        nearest_goal_id = knnsearch(centerLine', [curr_state(1), curr_state(3)]);
        lookahead_goal_id = nearest_goal_id + lookahead;
        if lookahead_goal_id > size(centerLine, 2)
            lookahead_goal_id = size(centerLine, 2);
        end
        disp("Tracking point (out of 246) and iteration: ");
        disp([nearest_goal_id, i]);
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
            [~,Y1]=ode45(@(t,x)bike(t,x,T,U_sim),T,curr_state);
            simulation_divergence = Y(end,:)-Y1(end,:);
        else
            [~,Y]=ode45(@(t,x)bike(t,x,T,U_sim),T,curr_state);
        end
        curr_state = Y(end,:);
        Y_checking = [Y_checking; [curr_state, goal_state, error, simulation_divergence]];
        prev_error = error;
        if passedGoal(curr_state(1), curr_state(3))
            break
        end
    end
end

function bool_out = passedGoal(x, y)
    bool_out = 0;
    if 61*x + 220*y + 269697 < 0
        bool_out = 1;
    end
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