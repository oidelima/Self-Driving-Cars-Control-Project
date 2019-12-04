function [L, dL_dU] = hilllossfun(curr_state, curr_input, hillMap) % TestTrack, cline_interp) %, goal_x, goal_y)
    x = curr_state(1);
    v = curr_state(2);
    y = curr_state(3);

    %L = hill(x,y, TestTrack, cline_interp);
    L = hill(x,y,v,hillMap);
    
    hill_step = 1; % .05;
    dL_dxy = gradientHill(x,y,v,hillMap, hill_step);

    % Step 1 forward integrate curr_input to get numerical gradient of
    % dx/dU and dy/dU with central difference
    step_size_u1 = .05; % delta goes from -.5 to .5
    step_size_u2 = 50; % Fx goes from -5000 to 5000
    num_simulated = 5; % 5 % >= 3 or the program crashes in ode45
    u1 = curr_input(1);
    u2 = curr_input(2);
    


    dt = .01; % * num_simulated;
    
    output_0    = curr_state;
    output_u1p1 = curr_state;
    output_u1m1 = curr_state;
    output_u2p1 = curr_state;
    output_u2m1 = curr_state;
    
    for i=1:num_simulated
        output_0    = output_0 + (bike(output_0',[u1 u2])*dt)';
        output_u1p1 = output_u1p1 + (bike(output_u1p1',[u1+step_size_u1 u2])*dt)';
        output_u1m1 = output_u1m1 + (bike(output_u1m1',[u1-step_size_u1 u2])*dt)';
        output_u2p1 = output_u2p1 + (bike(output_u2p1',[u1 u2+step_size_u2])*dt)';
        output_u2m1 = output_u2m1 + (bike(output_u2m1',[u1 u2-step_size_u2])*dt)';
    end
    
    x_0 = output_0(1,1);
    x_u1p1 = output_u1p1(1,1);
    x_u1m1 = output_u1m1(1,1);
    x_u2p1 = output_u2p1(1,1);
    x_u2m1 = output_u2m1(1,1); 
    y_0 = output_0(1,3);
    y_u1p1 = output_u1p1(1,3);
    y_u1m1 = output_u1m1(1,3);
    y_u2p1 = output_u2p1(1,3);
    y_u2m1 = output_u2m1(1,3);

    
    % dx_du1
    cd1_dx_du1 = (x_u1p1 - x_0) / step_size_u1;
    cd2_dx_du1 = (x_0 - x_u1m1) / step_size_u1;
    dx_du1 = (cd1_dx_du1 + cd2_dx_du1) / 2;
    
    % dx_du2
    cd1_dx_du2 = (x_u2p1 - x_0) / step_size_u2;
    cd2_dx_du2 = (x_0 - x_u2m1) / step_size_u2;
    dx_du2 = (cd1_dx_du2 + cd2_dx_du2) / 2;
    
    % dy_du1
    cd1_dy_du1 = (y_u1p1 - y_0) / step_size_u1;
    cd2_dy_du1 = (y_0 - y_u1m1) / step_size_u1;
    dy_du1 = (cd1_dy_du1 + cd2_dy_du1) / 2;
    
    % dy_du2
    cd1_dy_du2 = (y_u2p1 - y_0) / step_size_u2;
    cd2_dy_du2 = (y_0 - y_u2m1) / step_size_u2;
    dy_du2 = (cd1_dy_du2 + cd2_dy_du2) / 2;
    
    
    % dL_dU
    dL_dU = dL_dxy * [dx_du1 dx_du2; dy_du1 dy_du2];
    
end


function z = hill(x,y,v,hillMap)
    if y < -200 || y > 800 || x < 200 || x > 1500
        z = hill(1,1);
    else
        x = round(x);
        y = round(y);
        z = hillMap(y+201,x-199);
    end
    z = z/10000000000 + 1000*(10-v)^200; % 10000000000*(1/v)^4;
end





function dg = gradientHill(x,y,v,hillMap, step_size) % TestTrack,cline_interp,step_size)
    %f = @(x,y) hill(x,y,TestTrack,cline_interp) ;
    % f = @(x,y) hill2(x,y,TestTrack,cline_interp) ;
    f = @(x,y) hill(x,y,v,hillMap) ;
    [df_dx, df_dy] = brianNumGradient(f,x,y,step_size);
    dg = [df_dx df_dy];
end

function [df_dx, df_dy] = brianNumGradient(f,x,y,step_size)
    df_dx = (f(x+step_size,y) - f(x,y)) / step_size;
    df_dx = df_dx + ((f(x,y) - f(x-step_size,y)) / step_size);
    df_dx = df_dx / 2;
    
    df_dy = (f(x,y+step_size) - f(x,y)) / step_size;
    df_dy = df_dy + ((f(x,y) - f(x,y-step_size)) / step_size);
    df_dy = df_dy / 2;
end

function dzdt=bike(x,U)
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
delta_f = U(1,1);
F_x = U(1,2);

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



