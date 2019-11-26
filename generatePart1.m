clear
load TestTrack.mat
global n dt leftBound rightBound

% Required timestep
setGlobaldt(0.01);

% States: z = [x; u; y; v; psi; r]
z0 = [287; 5; -176; 0; 2; 0];

% Bounds on states and inputs
ub_0 = repmat([1500, 5000, 900, 5000, 3, 10], n, 1);
ub_1 = repmat([0.5, 5000], n-1, 1);
lb_0 = repmat([200, 0, -200, -5000, -3, -10], n, 1);
lb_1 = repmat([-0.5, -5000], n-1, 1);

ub = encodeColocationVector(ub_0(:,1), ub_0(:,2), ub_0(:,3), ub_0(:,4), ub_0(:,5), ub_0(:,6), ub_1(:,1), ub_1(:,2));
lb = encodeColocationVector(lb_0(:,1), lb_0(:,2), lb_0(:,3), lb_0(:,4), lb_0(:,5), lb_0(:,6), lb_1(:,1), lb_1(:,2));

% Track boundaries
leftBound = TestTrack.bl;
rightBound = TestTrack.br;

% Number of timesteps to use this iteration
setGlobaln(5000);

% Vector of inputs to be delivered for part 1; size (n, 2):
% [delta1, Fx1; delta2, Fx2; ...]
ROB535ControlsProjectpart1input = zeros(n,2);

% Fmincon
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
                       'SpecifyObjectiveGradient',true) ;
x0=zeros(1,8*n-2)';
cf=@costfun
nc=@nonlcon
z=fmincon(cf,x0,[],[],[],[],lb',ub',nc,options);

Y0=reshape(z(1:6*n),6,n)';
U=reshape(z(6*n+1:end),2,n-1);
u=@(t) [interp1(0:dt:(n-2)*dt,U(1,:),t,'previous','extrap');...
        interp1(0:dt:(n-2)*dt,U(2,:),t,'previous','extrap')];
[T1,Y1]=ode45(@(t,x) odefun(x,u(t)),[0:dt:(n-1)*dt],[0 0 0]);
[T2,Y2]=ode45(@(t,x) odefun(x,u(t)),[0:dt:(n-1)*dt],[0 0 -0.01]);
plot(Y0(:,1),Y0(:,2),Y1(:,1),Y1(:,2),Y2(:,1),Y2(:,2))
hold on
plot(0,0,'x');
legend('fmincon trajectory','ode45 trajectory using x0 = [0;0;0]',...
    'ode45 trajectory using x0 = [0;0;-0.01]','Buffered Obstacle','Start');
xlabel('x');
ylabel('y');

function [g,h,dg,dh]=nonlcon(z)
    global n dt leftBound rightBound
    [x, u, y, v, psi, r, delta_f, F_x] = decodeColocationVector(z);
    
    % Inequality Track Constraint, n x 1
    g = inequalityConstraintTrack([x,y]', leftBound, rightBound)';
    
    % Inequality Constraint Gradient, (8n - 2) x n
    dg_dxy = torGradient(@inequalityConstraintTrack, [x,y]', 0.05, leftBound, rightBound)';
    dg_c = mat2cell([dg_dxy(:,1), zeros(n,1), dg_dxy(:,2), zeros(n,3)]', [6], ones(1,n));
    dg = [blkdiag(dg_c{:}); zeros((n-1)*2, n)];
    
    % Equality Dynamics Constraint, (6n) x 1
    % Equality Constraint Gradient, (8n - 2) x (6n)
    % size of dh must be (8*n - 2), (6*n) = Transpose((no. of time steps * no. of states) x no. of elements in 'z') ;
    h = reshape(z(1:6),[6,1]);
    dh = zeros((8*n - 2), (6*n));
    dh(1:6,1:6) = eye(6);
    for i = 1:(n-1)
        curr_state = [x(i+1);u(i+1);y(i+1);v(i+1);psi(i+1);r(i+1)];
        prev_state = [x(i);u(i);y(i);v(i);psi(i);r(i)];
        prev_input = [delta_f(i);F_x(i)];
        h(6*i + 1:6*i + 6) = curr_state - prev_state - dt*bike_inst(prev_state, prev_input);
        dh(6*i + 1:6*i + 6, 6*i + 1:6*i + 6) = eye(6);
        dh(6*i - 5:6*i, 6*i + 1:6*i + 6) = -eye(6) - dt*jacobianDynamics(prev_state, prev_input, 1)';
        dh(6*n + 2*i - 1:6*n + 2*i, 6*i + 1:6*i + 6) = - dt*jacobianDynamics(prev_state, prev_input, 2)';
    end
    
end

function dzdt=bike_inst(x,u)
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

delta_f=u(1);
F_x=u(2);

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

function [J, dJ] = costfun(z)
    global n dt
    [x, u, y, v, psi, r, delta_f, F_x] = decodeColocationVector(z);
    input_cost = 0.25;
    J = norm(x - 1471.899)^2 + norm(y - 817.773)^2 + input_cost*(norm(delta_f)^2 + norm(F_x)^2);
    % size of dJ must be 1 x 603 (1 x no. of elements in 'z')
    zero_vec = zeros(size(x,1),1);
    dJ = encodeColocationVector(2*(x - 1471.899), zero_vec, 2*(y - 817.773), zero_vec, zero_vec, zero_vec, input_cost*2*delta_f, input_cost*2*F_x)';
end

% HW 3 Trajectory Generation:
%{
b = 1.5 ; 
L = 3 ;
dt=0.05;

%remember the format for z is as follows:
%z=[x0 y0 th0 x1 y1 th1 ... xn yn thn u0 d0 ... u(n-1) d(n-1)]';
nsteps = 121;
%3.1
ub = [];
ub2 = [];

lb = [];
lb2 = [];

for i = 0:120
    ub(3*i + 1) = 8;
    lb(3*i + 1) = -1;
    ub(3*i + 2) = 3;
    lb(3*i + 2) = -3;
    ub(3*i + 3) = pi/2;
    lb(3*i + 3) = -pi/2;
end
l1 = length(ub);
for j = 0:119
    ub2(2*j + 1) = 1;
    lb2(2*j + 1) = 0;
    ub2(2*j + 2) = 0.5;
    lb2(2*j + 2) = -0.5;
end
l2 = length(ub2);
ub = [ub,ub2]';
lb = [lb,lb2]';

%3.4
%%%%%%%%%%%%%%%% no need to change these lines  %%%%%%%%%%%%%%%%%%%%%%
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
                       'SpecifyObjectiveGradient',true) ;
x0=zeros(1,5*nsteps-2);
cf=@costfun
nc=@nonlcon
z=fmincon(cf,x0,[],[],[],[],lb',ub',nc,options);

Y0=reshape(z(1:3*nsteps),3,nsteps)';
U=reshape(z(3*nsteps+1:end),2,nsteps-1);
u=@(t) [interp1(0:dt:119*dt,U(1,:),t,'previous','extrap');...
        interp1(0:dt:119*dt,U(2,:),t,'previous','extrap')];
[T1,Y1]=ode45(@(t,x) odefun(x,u(t)),[0:dt:120*dt],[0 0 0]);
[T2,Y2]=ode45(@(t,x) odefun(x,u(t)),[0:dt:120*dt],[0 0 -0.01]);
plot(Y0(:,1),Y0(:,2),Y1(:,1),Y1(:,2),Y2(:,1),Y2(:,2))
theta = 0:0.01:2*pi;
hold on
plot((0.7*cos(theta)+3.5),(0.7*sin(theta)-0.5))
hold on
plot(0,0,'x');
legend('fmincon trajectory','ode45 trajectory using x0 = [0;0;0]',...
    'ode45 trajectory using x0 = [0;0;-0.01]','Buffered Obstacle','Start');
ylim([-2,2]);
xlim([-1,8]);
xlabel('x');
ylabel('y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3.2
%z=[x0 y0 th0 x1 y1 th1 ... xn yn thn u0 d0 ... u(n-1) d(n-1)]';
function [g,h,dg,dh]=nonlcon(z)
    b = 1.5 ; 
    L = 3 ;
    dt=0.05;
    % size of g must be 121 x 1 (no.of time steps);
    l1 = 363;
    l2 = 240;
    [x,y,phi,u,del] = decode(z);
    g = (0.7)^2 - (x - 3.5).^2 - (y - (-0.5)).^2;
    dg_dx = 7 - 2*x;
    dg_dy = -2*y - 1;
    size(dg_dx)
    size(dg_dy)
    %dg = [repmat(encode(dg_dx, dg_dy, zeros(121,1), zeros(120,1), zeros(120,1)), 1, 121)];
    dg_c = mat2cell([dg_dx, dg_dy, zeros(121,1)]', [3], ones(1,121))
    dg_c{1}
    dg_c{5}
    dg = [blkdiag(dg_c{:}); zeros(240, 121)];
    
    % size of dg must be 603 x 121 = Transpose(no. of time steps x no. of elements in 'z');
    % size of h must be 363  1 ((no. of time steps * no. of states) x 1)
    % size of dh must be 603 x 363 = Transpose((no. of time steps * no. of states) x no. of elements in 'z') ;
    h = reshape(z(1:3),[3,1]);
    dh = zeros(603, 363);
    dh(1:3,1:3) = eye(3);
    size(x)
    size(y)
    size(phi)
    for i = 1:120
        h = [h;[x(i+1);y(i+1);phi(i+1)] - [x(i);y(i);phi(i)] - 0.05*odefun([x(i);y(i);phi(i)],[u(i);del(i)])];
        dh(3*i + 1:3*i + 3, 3*i + 1:3*i + 3) = eye(3);
        dh(3*i - 2:3*i, 3*i + 1:3*i + 3) = -eye(3) - 0.05 * [0, 0, -u(i)*sin(phi(i))-(b/L)*u(i)*tan(del(i))*cos(phi(i));0, 0, u(i)*cos(phi(i))-(b/L)*u(i)*tan(del(i))*sin(phi(i)); 0, 0, 0]';
        dh(363 + 2*i - 1:363 + 2*i, 3*i + 1:3*i + 3) = - 0.05 * [cos(phi(i))-(b/L)*tan(del(i))*sin(phi(i)), -(b/L)*u(i)*(sec(del(i))^2)*sin(phi(i)); sin(phi(i))+(b/L)*tan(del(i))*cos(phi(i)), (b/L)*u(i)*(sec(del(i))^2)*cos(phi(i)); (1/L)*tan(del(i)), (u(i)/L)*sec(del(i))^2]';
    end
    
end

%3.3
function [J, dJ] = costfun(z)
    % size of J must be 1 x 1
    size(z)
    [x,y,phi,u,del] = decode(z);
    J = norm(x - 7)^2 + norm(y)^2 + norm(phi)^2 + norm(u)^2 + norm(del)^2;
    % size of dJ must be 1 x 603 (1 x no. of elements in 'z')
    dJ = encode(2*(x-7), 2*y, 2*phi, 2*u, 2*del)';
end


function [dx] = odefun(x,u)
    b = 1.5 ; 
    L = 3 ;
    dx = [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
          u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
          u(1)*tan(u(2))/L] ;
end

function [x,y,psi,u,del] = decode(z)
    l1 = 363;
    l2 = 240;
    z = reshape(z, [l1+l2,1]);
    x = z(1:3:l1);
    y = z(2:3:l1);
    psi = z(3:3:l1);
    u = z(l1+1:2:l1+l2);
    del = z(l1+2:2:l1+l2);
end

function z = encode(x,y,psi,u,del)
    l1 = 363;
    l2 = 240;
    z1 = reshape([x,y,psi]',[l1,1]);
    z2 = reshape([u,del]',[l2,1]);
    z = [z1;z2];
end
%}