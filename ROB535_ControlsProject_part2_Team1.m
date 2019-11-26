%function [sol_2, FLAG_terminate] = ROB535_ControlsProject_part2_Team1 (TestTrack,Xobs_seen,curr_state)


if curr_state == 
    FLAG_terminate = 1
else
    FLAG_terminate = 0
end


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
xsize = 6*51;
zsize = Ndec;

% Staates : z  &  Input : u 
%z = zeros(6*51 + 2*5, 1 );



% State Constraints 

[leftlimit, rightlimit] = bound(npred,states,leftbound, rightbound);
state_lb = -inf(Ndec,1);
state_ub = inf(Ndec,1);

for i = 1 : npred
    state_lb(6* (i-1) +1) =  min(leftlimit(1,i), rightlimit(1,i));
    state_lb(6* (i-1) +3) =  min(leftlimit(2,i), rightlimit(2,i));
    
    state_ub(6* (i-1) +1) = max(leftlimit(1,i), rightlimit(1,i));
    state_ub(6* (i-1) +3) = max(leftlimit(2,i), rightlimit(2,i));
end




% Input constraints 
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
% !! have to call the result of part 1 for calling Y_ref / U_ref !! 
    STATE_tilde = urr_state - Y_ref;
    [Aeq,beq] = eq_cons(npred, xsize, zsize, A,B, STATE_tilde);
    [Lb,Ub] = bound_cons(U_ref);
    
    
    sol = quadprog(H,zeros(53,1),[],[],Aeq,beq,Lb,Ub);


    %extract control input from solution   
    %forward integrate using nonlinear model, STATE(i), the extracted control input and ode45
    uin=sol(xsize+1 : xsize+2,1) + U_ref;



%% FLAG terminate is a binary flag set by the function to indicate when to stop the simulation. Setting this flag to 1 at a planning iteration indicates that, after integrating forward vehicle dynamics using the control inputs from sol 2 for 0.5 second, the vehicle would have reached the goal. If this flag is set to 0, then it implies that the car would not reach the goal after forward integrating for 0.5 second. Note that, in the event that FLAG terminate is never set to 1, the simulation would stop after 20 minutes so that we can process all the teams files in a timely manner.





function [Aeq,Beq] = eq_cons(npred, xsize, zsize, A,B, curr_state)

Aeq=zeros(xsize,zsize);
Beq=zeros(xsize,1);

for i = 1 : npred
    Am = A(i);
    Bm = B(i);
    
    Aeq(6*(i)+1 :6*(i)+6 , 6*(i-1)+1:6*(i-1)+6) = Am ;
    Aeq(6*(i)+1: 6*(i)+6, 6*(i)+1 :6*(i)+6 ) = Aeq(6*(i)+1: 6*(i)+6, 6*(i)+1 :6*(i)+6) -eye(6);  
    
    Aeq(6*(i)+1:6*(i)+6 ,xsize+2*(i-1)+1:xsize+2*(i-1)+2) = Bm;
end

Aeq(1:6,1:6) = eye(6);
Beq(1:6) = reshape(curr_state,6,1);

end



function[leftlimit, rightlimit] = bound(npred,states,leftbound, rightbound)
    
    % Output : leftlimit and the right limit of the road for the current
    % location
     
    xypos = zeros(npred+1,2);
    leftlimit = -inf(size(states))
    rightlimit = inf(size(states))
    
    
    for i = 1 : npred+1
        xypos(i,1) = states(6 * (i-1) +1);
        xypos(i,2) = states(6 * (i-1) +3);
        Idleft_1 = knnsearch(leftBound', xypos(i,:));
        Idright_1 = knnsearch(rightBound', xypos(i,:));
        
        % find a left and right pair of points to do a linear interpolation
        % in order to get the exact leftbound and rightbound of a given
        % position. Linear interpolation is done by using a "dot product".
        
        if dot(xypos(i,:) - leftbound(:,Idleft_1),leftbound(:, Idleft_1 +1)-leftbound(:, Idleft_1)) >= 0 
            leftlimit = leftbound(:,Idleft_1) + dot(xypos(i,:) - leftbound(:,Idleft_1),leftbound(:, Idleft_1 +1)-leftbound(:, Idleft_1));
        else
            leftlimit = leftbound(:,Idleft_1) + dot(xypos(i,:) - leftbound(:,Idleft_1),leftbound(:, Idleft_1-1)-leftbound(:, Idleft_1));
        end
        
        if dot(xypos(i,:) - rightbound(:,Idright_1),rightbound(:,Idright_1 +1)-xypos(i,:)) >= 0 
            rightlimit = dot(xypos(i,:) - rightbound(:,Idright_1),rightbound(:,Idright_1 +1)-xypos(i,:));
        else
            rightlimit = dot(xypos(i,:) - rightbound(:,Idright_1),rightbound(:,Idright_1 -1)-xypos(i,:));
        end
    end 
 
end


%function restricting ith x and y with 1st order function
% for a and b 

function[a, b] = first_ineq(zsize,npred,states,leftbound,rightbound)
   
    a = zeros(2*(npred+1), zsize);
    b = zeros(2*(npred+1),1);

     
     for i = 1 : npred+1
        xpos = states(6 * (i-1) +1);
        ypos = states(6 * (i-1) +3);
        Idleft_1 = knnsearch(leftBound', [xpos,ypos]');
        Idright_1 = knnsearch(rightBound',[xpos,ypos]');
        
        left_1 = leftbound(:,Idleft_1);
        right_1 = rightbound(:,Idright_1);
        
        if dot([xpos,ypos]' - left_1,leftbound(:, Idleft_1 +1)-left_1) >= 0 
            left_2 = left_1  + dot([xpos,ypos]' - lleft_1 ,leftbound(:, Idleft_1 +1)-left_1 );
        else
            left_2 = left_1  + dot([xpos,ypos]' - left_1 ,leftbound(:, Idleft_1-1)-left_1 );
        end
        
        if dot([xpos,ypos]' - left_1,rightbound(:,Idright_1 +1)-left_1) >= 0 
            right_2 = right_1+ dot([xpos,ypos]' - right_1,rightbound(:,Idright_1 +1)-right_1);
        else
            right_2 = right_1+ dot([xpos,ypos]' - right_1,rightbound(:,Idright_1 -1)-right_1);
        end
        
        left_slope = (left_1(1) - left_2(1)) / (left_1(2) - left_2(2));
        right_slope = (right_1(1) - right_2(1)) / (right_1(2) - right_2(2));
        
        a(2*(i-1) +1, 6 * (i-1) +1) = - left_slope;
        a(2*(i-1) +1, 6 * (i-1) +3) = 1;
        a(2*(i-1) +2, 6 * (i-1) +1) = - left_slope;
        a(2*(i-1) +1, 6 * (i-1) +3) = -1;
        
        b(2*(i-1) +1) = left_slope * (-left_1(1)) + left_1(2);
        b(2*(i-1) +2) = right_slope * (left_1(1)) - left_1(2);
        
     end       
 end
     

% avoiding obstacles by elongating(adding) some rows a and b constraints % 

function avoiding_obstacles(Xobs_seen)
    obs_num = length(Xobs_seen)

end
