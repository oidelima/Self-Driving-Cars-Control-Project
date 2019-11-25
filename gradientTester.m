clear
format long
load TestTrack.mat
curr_state =[287,5,-176,0,2,0];

% Constants
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

dt = .01;

% Symbollic state
syms x_sym u_sym y_sym v_sym psi_sym r_sym

% Symbollic input
syms delta_f_sym F_x_sym

% Symbollic goal states
syms x_goal_sym y_goal_sym

%slip angle functions in degrees
a_f=rad2deg(delta_f_sym-atan2(v_sym+a*r_sym,u_sym));
a_r=rad2deg(-atan2(v_sym-b*r_sym,u_sym));

%Nonlinear Tire Dynamics
phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

F_zf=b/(a+b)*m*g;
F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

F_zr=a/(a+b)*m*g;
F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

F_total=simplify(sqrt((Nw*F_x_sym)^2+(F_yr^2))); % need this for if
F_max=0.7*m*g; % need this for if

% Can still deal with max, not sure if necessary
% %if F_total>F_max
%     F_x=F_max/F_total*F_x_sym;
%   
%     F_yr=F_max/F_total*F_yr;
% end

%vehicle dynamics
dz_dt= [u_sym*cos(psi_sym)-v_sym*sin(psi_sym);...
          (-f*m*g+Nw*F_x_sym-F_yf*sin(delta_f_sym))/m+v_sym*r_sym;...
          u_sym*sin(psi_sym)+v_sym*cos(psi_sym);...
          (F_yf*cos(delta_f_sym)+F_yr)/m-u_sym*r_sym;...
          r_sym;...
          (F_yf*a*cos(delta_f_sym)-F_yr*b)/Iz]



dx_dt = dz_dt(1)
du_dt = dz_dt(2)
dy_dt = dz_dt(3)
dv_dt = dz_dt(4)
dpsi_dt = dz_dt(5)
dr_dt = dz_dt(6)

% Calc next x and y as function of previous symbollic 
x = x_sym + dx_dt * dt
y = y_sym + dy_dt * dt

% u = u_sym + du_dt * dt
% r = r_sym + dr_dt * dt
% psi = psi_sym + dpsi_dt * dt
% v = v_sym + dv_dt
% 
% du_dU = simplify(jacobian(u,[delta_f_sym F_x_sym]))
% dr_dU = simplify(jacobian(r,[delta_f_sym F_x_sym]))
% dpsi_dU = simplify(jacobian(psi,[delta_f_sym F_x_sym]))
% dv_dU = simplify(jacobian(v,[delta_f_sym F_x_sym]))

dx_dU = simplify(jacobian(x,[delta_f_sym F_x_sym]))
dy_dU = simplify(jacobian(y,[delta_f_sym F_x_sym]))

%%

%LossFun = @(curr_state, x_goal, y_goal)  % include v^2?
LossFun_sym = simplify((x - x_goal_sym)^2 + (y - y_goal_sym)^2)

dL_dU_sym = simplify(jacobian(LossFun_sym, [delta_f_sym F_x_sym]))
%dL_dU = subs(dL_dU_sym, [x_sym u_sym y_sym v_sym psi_sym r_sym], curr_state)





