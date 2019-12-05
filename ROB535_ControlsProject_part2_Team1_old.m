function [sol_2, FLAG_terminate] = ROB535_ControlsProject_part2_Team1_old (TestTrack,Xobs_seen,curr_state)
    leftBound = TestTrack.bl;
    rightBound = TestTrack.br;
    centerLine = TestTrack.cline;
    centerLineTheta = TestTrack.theta;
    FLAG_terminate = 0;
    
    % Times in seconds
    total_time = .7;
    dt = .01;
    
    kp1 = 0.5;
    kp2 = 1000;
    
    speed = 6;
    
    % Get cnormals
    n_elements = size(centerLine, 2);
    projectors = (centerLine(:,2:n_elements) - centerLine(:,1:n_elements-1));
    projectors = projectors ./ vecnorm(projectors);
    Cnormals = [0, 1; -1, 0] * projectors;
    Cnormals(:,n_elements) = Cnormals(:,n_elements - 1);
    
    nsteps = total_time/dt;
    
    U = zeros(nsteps,2);
    %U(1,:) = [0 0]
    
    for i=2:nsteps+1
        nearest_goal_id = knnsearch(centerLine', [curr_state(1), curr_state(3)]);% + lookahead;
        %disp("Tracking point (out of 246) and iteration: ");
        disp([nearest_goal_id, i]);
        nearest_goal = centerLine(:, nearest_goal_id);
        goal_state = [nearest_goal(1), speed, nearest_goal(2), 0, 0, 0];
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
        %U = [U;[delta_f, Fx]];
        U(i,:) = [delta_f, Fx];
        
        U_sim = U(i-1:i, :);
        x0 = curr_state;
        [~,Y]=ode45(@(t,x)bike(t,x,T,U_sim),T,x0);
        curr_state = Y(end,:);
        
        if close_to_end(curr_state,TestTrack)
            FLAG_terminate = 1;
        end
        
        %Y_checking = [Y_checking; [curr_state, goal_state, error, ]];
    end
    
    sol_2 = U(2:end);

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
    is_close = ( dist < 1 );
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
% % (as an N × 2 vector where the first column is ? and the second column is Fx
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
