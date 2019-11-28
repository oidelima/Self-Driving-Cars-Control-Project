function [Y_checking, U] = PDControlPart1(nsteps, x0, gains, speed)
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
    
    kp1 = 0.5;
    kp2 = 1000;
    
    % Required timestep
    setGlobaldt(0.01);

    curr_state = x0;
    n_elements = size(centerLine, 2);
    projectors = (centerLine(:,2:n_elements) - centerLine(:,1:n_elements-1));
    projectors = projectors ./ vecnorm(projectors);
    Cnormals = [0, 1; -1, 0] * projectors;
    Cnormals(:,n_elements) = Cnormals(:,n_elements - 1);

    U = [0, 0];
    Y_checking = [];

    for i=1:nsteps
        nearest_goal_id = knnsearch(centerLine', [curr_state(1), curr_state(3)]);% + lookahead;
        disp("Tracking point (out of 246) and iteration: ");
        disp([nearest_goal_id, i]);
        nearest_goal = centerLine(:, nearest_goal_id);
        goal_state = [nearest_goal(1), 5, nearest_goal(2), 0, 0, 0];
        error = goal_state - curr_state;
        goal_state(5) = centerLineTheta(:, nearest_goal_id) -.1*dot(Cnormals(:, nearest_goal_id), [error(1);error(3)], 1); %atan2(error(3),error(1));
        error = goal_state - curr_state;
        delta_f = kp1 * error(5);
        if delta_f > 0.5
            delta_f = 0.5;
        end
        if delta_f < -0.5
            delta_f = -0.5;
        end
        Fx = kp2 * error(2);
        %Solve for trajectory
        T=(i-1)*0.01:0.01:i*0.01;
        U = [U;[delta_f, Fx]];
        U_sim = U(end-1:end, :);
        x0 = curr_state;
        [~,Y]=ode45(@(t,x)bike(t,x,T,U_sim),T,x0);
        curr_state = Y(end,:);
        Y_checking = [Y_checking; [curr_state, goal_state, error, ]];
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