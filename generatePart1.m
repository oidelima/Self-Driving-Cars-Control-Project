clear
load TestTrack.mat
global n dt leftBound rightBound F_max

% Number of timesteps to use this iteration
setGlobaln(1000);

% Required timestep
setGlobaldt(0.01);

% Store precomputed symbolic jacobian saturation threshold
m = 1400;
g = 9.806;
F_max = 0.7*m*g;

% States: z = [x; u; y; v; psi; r]
z0 = repmat([287, 5, -176, 0, 2, 0], n, 1);
u0 = repmat([0 0], n-1, 1);
z0_colocated = encodeColocationVector(z0(:,1), z0(:,2), z0(:,3), z0(:,4), z0(:,5), z0(:,6), u0(:,1), u0(:,2));

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

% Vector of inputs to be delivered for part 1; size (n, 2):
% [delta1, Fx1; delta2, Fx2; ...]
ROB535ControlsProjectpart1input = zeros(n,2);

% Fmincon
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
                       'SpecifyObjectiveGradient',true,...
                       'MaxIterations',1500,...
                       'Display', 'iter') ;
x0=z0_colocated
cf=@costfun
nc=@nonlcon
z=fmincon(cf,x0,[],[],[],[],lb',ub',nc,options);

[x, u, y, v, psi, r, delta, Fx] = decodeColocationVector(z)
out_input = [delta,Fx];
legend('fmincon trajectory','ode45 trajectory using x0 = [0;0;0]',...
    'ode45 trajectory using x0 = [0;0;-0.01]','Buffered Obstacle','Start');
xlabel('x');
ylabel('y');

function [g,h,dg,dh]=nonlcon(z)
    global n dt leftBound rightBound F_max
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
    h = reshape(z(1:6),[6,1]) - [287, 5, -176, 0, 2, 0]';
    dh = zeros((8*n - 2), (6*n));
    dh(1:6,1:6) = eye(6);
    for i = 1:(n-1)
        curr_state = [x(i+1);u(i+1);y(i+1);v(i+1);psi(i+1);r(i+1)];
        prev_state = [x(i);u(i);y(i);v(i);psi(i);r(i)];
        prev_input = [delta_f(i);F_x(i)];
        h(6*i + 1:6*i + 6) = curr_state - prev_state - dt*bikeInst(prev_state, prev_input);
        dh(6*i + 1:6*i + 6, 6*i + 1:6*i + 6) = eye(6);
        if F_total(prev_state,prev_input) > F_max
            J1used = @J1sat;
            J2used = @J2sat;
        else
            J1used = @J1;
            J2used = @J2;
        end
        dh(6*i - 5:6*i, 6*i + 1:6*i + 6) = -eye(6) - dt*J1used(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_input(1), prev_input(2))';
        dh(6*n + 2*i - 1:6*n + 2*i, 6*i + 1:6*i + 6) = - dt*J2used(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_input(1), prev_input(2))';
    end
    
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

function F = F_total(x,u)
    % States: z = [x; u; y; v; psi; r]
    F = (4*u(2)^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*x(4) - 1.4500*x(6), x(2)) + 1.6000*atan(15.4699*atan2(1*x(4) - 1.4500*x(6), x(2)))))^2)^(1/2);
end

function J = J1(x,u,y,v,psi_,r,delta_f,F_x)
    J = [0,cos(psi_), 0,-sin(psi_), - v*cos(psi_) - u*sin(psi_),0;
        0,-(4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f)*((40.2216*(1.3500*r + v))/((1.3500*r + v)^2 + u^2) - (24.7518*(1.3500*r + v))/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1), 0,r + (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f)*((40.2216*u)/((1.3500*r + v)^2 + u^2) - (24.7518*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1),0,v + (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f)*((54.2992*u)/((1.3500*r + v)^2 + u^2) - (33.4149*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1);
        0,sin(psi_),0,cos(psi_),u*cos(psi_) - v*sin(psi_),0;
        0, - r - (3.9714*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*(1.4500*r - v))/((1.4500*r - v)^2 + u^2) - (24.7518*(1.4500*r - v))/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1) + (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((40.2216*(1.3500*r + v))/((1.3500*r + v)^2 + u^2) - (24.7518*(1.3500*r + v))/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1), 0, - (3.9714*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*u)/((1.4500*r - v)^2 + u^2) - (24.7518*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1) - (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((40.2216*u)/((1.3500*r + v)^2 + u^2) - (24.7518*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1),0, - u + (3.9714*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((58.3214*u)/((1.4500*r - v)^2 + u^2) - (35.8901*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1) - (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((54.2992*u)/((1.3500*r + v)^2 + u^2) - (33.4149*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1);
        0,0, 0,0,0,1;
        0,(3.0229*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*(1.4500*r - v))/((1.4500*r - v)^2 + u^2) - (24.7518*(1.4500*r - v))/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1) + (3.0229*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((40.2216*(1.3500*r + v))/((1.3500*r + v)^2 + u^2) - (24.7518*(1.3500*r + v))/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1), 0,   (3.0229*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*u)/((1.4500*r - v)^2 + u^2) - (24.7518*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1) - (3.0229*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((40.2216*u)/((1.3500*r + v)^2 + u^2) - (24.7518*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1),0,- (3.0229*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((58.3214*u)/((1.4500*r - v)^2 + u^2) - (35.8901*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1) - (3.0229*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((54.2992*u)/((1.3500*r + v)^2 + u^2) - (33.4149*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1)];
end

function J = J1sat(x,u,y,v,psi_,r,delta_f,F_x)
    J = [0,cos(psi_), 0,-sin(psi_), - v*cos(psi_) - u*sin(psi_),0;
        0,- (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f)*((40.2216*(1.3500*r + v))/((1.3500*r + v)^2 + u^2) - (24.7518*(1.3500*r + v))/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1) + (3.5366e+08*F_x*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*(1.4500*r - v))/((1.4500*r - v)^2 + u^2) - (24.7518*(1.4500*r - v))/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000), 0,r + (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f)*((40.2216*u)/((1.3500*r + v)^2 + u^2) - (24.7518*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1) + (3.5366e+08*F_x*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*u)/((1.4500*r - v)^2 + u^2) - (24.7518*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000),0,v + (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f)*((54.2992*u)/((1.3500*r + v)^2 + u^2) - (33.4149*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1) - (3.5366e+08*F_x*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((58.3214*u)/((1.4500*r - v)^2 + u^2) - (35.8901*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000);
        0,sin(psi_), 0,cos(psi_),   u*cos(psi_) - v*sin(psi_),0;
        0, - r - (3.8165e+04*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*(1.4500*r - v))/((1.4500*r - v)^2 + u^2) - (24.7518*(1.4500*r - v))/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^0.5000) + (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((40.2216*(1.3500*r + v))/((1.3500*r + v)^2 + u^2) - (24.7518*(1.3500*r + v))/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1) + (8.1932e+11*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2*((40.2216*(1.4500*r - v))/((1.4500*r - v)^2 + u^2) - (24.7518*(1.4500*r - v))/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000), 0, - (3.8165e+04*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*u)/((1.4500*r - v)^2 + u^2) - (24.7518*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^0.5000) - (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((40.2216*u)/((1.3500*r + v)^2 + u^2) - (24.7518*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1) + (8.1932e+11*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2*((40.2216*u)/((1.4500*r - v)^2 + u^2) - (24.7518*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000),0, - u + (3.8165e+04*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((58.3214*u)/((1.4500*r - v)^2 + u^2) - (35.8901*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^0.5000) - (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((54.2992*u)/((1.3500*r + v)^2 + u^2) - (33.4149*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1) - (8.1932e+11*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2*((58.3214*u)/((1.4500*r - v)^2 + u^2) - (35.8901*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000);
        0,0, 0,0,0,1;
        0,(2.9049e+04*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*(1.4500*r - v))/((1.4500*r - v)^2 + u^2) - (24.7518*(1.4500*r - v))/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^0.5000) + (3.0229*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((40.2216*(1.3500*r + v))/((1.3500*r + v)^2 + u^2) - (24.7518*(1.3500*r + v))/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1) - (6.2363e+11*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2*((40.2216*(1.4500*r - v))/((1.4500*r - v)^2 + u^2) - (24.7518*(1.4500*r - v))/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000), 0,   (2.9049e+04*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((40.2216*u)/((1.4500*r - v)^2 + u^2) - (24.7518*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^0.5000) - (3.0229*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((40.2216*u)/((1.3500*r + v)^2 + u^2) - (24.7518*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1) - (6.2363e+11*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2*((40.2216*u)/((1.4500*r - v)^2 + u^2) - (24.7518*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000),0,- (2.9049e+04*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*((58.3214*u)/((1.4500*r - v)^2 + u^2) - (35.8901*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^0.5000) - (3.0229*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*((54.2992*u)/((1.3500*r + v)^2 + u^2) - (33.4149*u)/(((1.3500*r + v)^2 + u^2)*((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1))))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1) + (6.2363e+11*cos(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2*((58.3214*u)/((1.4500*r - v)^2 + u^2) - (35.8901*u)/((239.3166*atan2(1*v - 1.4500*r, u)^2 + 1)*((1.4500*r - v)^2 + u^2))))/(((40.2216*atan2(1*v - 1.4500*r, u) - 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))^2 + 1)*(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000)];
 
end

function J = J2(x,u,y,v,psi_,r,delta_f,F_x)
    J = [0,0;
        3.5547*sin(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f) + (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f)*(24.7518/((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1) - 40.2216))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1),0.0014;
        0,0;
        3.5547*sin(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f) - (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*(24.7518/((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1) - 40.2216))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1),0;
        0,0;
        2.5191*sin(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f) - (3.0229*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*(24.7518/((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1) - 40.2216))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1),0];
 
end

function J = J2sat(x,u,y,v,psi_,r,delta_f,F_x)
    J = [0,0;
        3.5547*sin(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f) + (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f)*(24.7518/((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1) - 40.2216))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1), 13.7284/(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^0.5000 - (54.9136*F_x^2)/(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000;
        0,0;
        3.5547*sin(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f) - (4.2656*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*(24.7518/((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1) - 40.2216))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1),                                         -(1.2722e+05*F_x*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))))/(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000;
        0,0;
        2.5191*sin(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*sin(delta_f) - (3.0229*cos(1.2000*atan(- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u)))*cos(delta_f)*(24.7518/((15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u))^2 + 1) - 40.2216))/((- 40.2216*delta_f + 1.6000*atan(15.4699*delta_f - 15.4699*atan2(1.3500*r + 1*v, u)) + 40.2216*atan2(1.3500*r + 1*v, u))^2 + 1),                                          (9.6831e+04*F_x*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u)))))/(4*F_x^2 + 2.1468e+07*sin(1.2000*atan(- 40.2216*atan2(1*v - 1.4500*r, u) + 1.6000*atan(15.4699*atan2(1*v - 1.4500*r, u))))^2)^1.5000];
end
%}

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