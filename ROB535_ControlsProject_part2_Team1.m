%function [sol_2, FLAG_terminate] = ROB535_ControlsProject_part2_Team1 (TestTrack,Xobs_seen,curr_state)

%%% where sol 1 is the vehicle?s trajectory.

%% sol 2 is a vector of control inputs that will be passed to forwardIntegrateControlInput. 
%As mentioned later in 6.2, sol 2 must have enough control inputs with time step 0.01 sec- onds to forward integrate vehicle dynamics for the next 0.5 second.

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


% (as an N × 2 vector where the first column is ? and the second column is Fx

%% linearize 

syms x1 x2 x3 x4 x5 x6 u1 u2
 
 
a_f=rad2deg(u1-atan2(x4+a*x6,x2));
a_r=rad2deg(-atan2((x4-b*x6),x2));
 
 
F_zf=b/(a+b)*m*g;
F_yf=F_zf*Dy*sin(Cy*atan(By*((1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy)))))+Svy;
 
F_zr=a/(a+b)*m*g;
F_yr=F_zr*Dy*sin(Cy*atan(By*((1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy)))))+Svy;
 
 
dzdx= [x2*cos(x5)-x4*sin(x5);...
          (-f*m*g+Nw*u2-F_yf*sin(u1))/m+x4*x6;...
          x2*sin(x5)+x4*cos(x5);...
          (F_yf*cos(u1)+F_yr)/m-x2*x6;...
          x6;...
          (F_yf*a*cos(u1)-F_yr*b)/Iz];
 
A = jacobian(dzdx,[x1;x2;x3;x4;x5;x6]);
B = jacobian(dzdx,[u1;u2]);
 
%% Decision Variables 
% ith # of state : 6, input : 2 ;
% horizon of 0.5 s & 0.01 s = 50 timestep
Ndec = 6*51 + 2*5;
npred = 50;

% Staates : z  &  Input : u 
%z = zeros(6*51 + 2*5, 1 );
%dzdt=bike(t,x,T,U)

%% Equality Constraints
% quadprog
% Euler integration  
xsize = 6*51;
zsize = Ndec;
Aeq=zeros(xsize,zsize);
Beq=zeros(xsize,1);


for i = 1 : npred
    Am = A(i);
    Bm = B(i);
    
    Aeq(6*(i)+1 :6*(i)+6 , 6*(i-1)+1:6*(i-1)+6) = Am ;
    Aeq(6*(i)+1:6*(i)+6, 6*(i-1)+1+6) = -ones(6,1);  
    
    Aeq(6*(i)+1:6*(i)+6 ,xsize+2*(i-1)+1:xsize+2*(i-1)+2) = Bm;
end

Aeq(1:6,1:6) = eye(6);
Beq(1:6) = reshape(curr_state,6,1);


%%Boundary Constraints
% 4-by-2 matrix 
    Lb = -inf(Ndec,1);
    Ub = inf(Ndec,1);
    j=0;
    for k = 1: npred
        Lb(xsize+2*(k-1) + 1) = 0 - U(1,i+j);
        Lb(xsize+2*(k-1) + 2) = -0.5 - U(2,i+j);
 
        Ub(xsize+2*(k-1) + 1) = 1 -U(1,i+j);
        Ub(xsize+2*(k-1) + 2) = 0.5 -U(2,i+j);
        j = j +1;
    end


%% quadatic penalty




%% FLAG terminate is a binary flag set by the function to indicate when to stop the simulation. Setting this flag to 1 at a planning iteration indicates that, after integrating forward vehicle dynamics using the control inputs from sol 2 for 0.5 second, the vehicle would have reached the goal. If this flag is set to 0, then it implies that the car would not reach the goal after forward integrating for 0.5 second. Note that, in the event that FLAG terminate is never set to 1, the simulation would stop after 20 minutes so that we can process all the teams files in a timely manner.










%
