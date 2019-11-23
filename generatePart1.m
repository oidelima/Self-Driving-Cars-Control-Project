clear
load TestTrack.mat
load gridInequal.mat

% Required timestep
dt = 0.01;

% States: z = [x; u; y; v; psi; r]
z0 = [287; 5; -176; 0; 2; 0];

% Number of timesteps to use this iteration
n = 500;

% Vector of inputs to be delivered for part 1; size (n, 2):
% [delta1, Fx1; delta2, Fx2; ...]
ROB535ControlsProjectpart1input = zeros(n,2);

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