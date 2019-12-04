function [L, dL_dU] = lossfun(curr_state, curr_input, TestTrack, cline_interp) %, goal_x, goal_y)
    x = curr_state(1);
    y = curr_state(3);

    %L = hill(x,y, TestTrack, cline_interp);
    L = hill2(x,y, TestTrack, cline_interp);
    
    hill_step = .1; % .05;
    dL_dxy = gradientHill(x,y,TestTrack,cline_interp, hill_step);

    % Step 1 forward integrate curr_input to get numerical gradient of
    % dx/dU and dy/dU with central difference
    step_size_u1 = .05; % delta goes from -.5 to .5
    step_size_u2 = 50; % Fx goes from -5000 to 5000
    num_simulated = 3; % 5 % >= 3 or the program crashes in ode45
    u1 = curr_input(1);
    u2 = curr_input(2);
    
    % Simulate the resulting Y from the inputs
    [output_0,~] = forwardIntegrateControlInput([u1*ones(num_simulated,1), u2*ones(num_simulated,1)], curr_state);
    
    % Shift u1 by step size and forward integrate
    [output_u1p1,~] = forwardIntegrateControlInput([(u1+step_size_u1)*ones(num_simulated,1), u2*ones(num_simulated,1)], curr_state);
    [output_u1m1,~] = forwardIntegrateControlInput([(u1-step_size_u1)*ones(num_simulated,1), u2*ones(num_simulated,1)], curr_state);
    
    % Shift u2 by step size and forward integrate
    [output_u2p1,~] = forwardIntegrateControlInput([u1*ones(num_simulated,1), (u2+step_size_u2)*ones(num_simulated,1)], curr_state);
    [output_u2m1,~] = forwardIntegrateControlInput([u1*ones(num_simulated,1), (u2-step_size_u2)*ones(num_simulated,1)], curr_state);
    
    % Assign x's and y's from the last step of the simulation
    x_0 = output_0(end,1);
    x_u1p1 = output_u1p1(end,1);
    x_u1m1 = output_u1m1(end,1);
    x_u2p1 = output_u2p1(end,1);
    x_u2m1 = output_u2m1(end,1); 
    y_0 = output_0(end,3);
    y_u1p1 = output_u1p1(end,3);
    y_u1m1 = output_u1m1(end,3);
    y_u2p1 = output_u2p1(end,3);
    y_u2m1 = output_u2m1(end,3); 
    
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

function z = hill2(x,y, TestTrack, cline_interp)
%     pos = [x;y];
%     leftBound = TestTrack.bl;
%     rightBound = TestTrack.br;
%     [leftNormals, rightNormals] = calcNormals(leftBound, rightBound);
%     Idleft = knnsearch(leftBound', pos');
% 	Idright = knnsearch(rightBound', pos');
% 
%     rightCorr = rightBound(:, Idright');
%     leftCorr = leftBound(:, Idleft');
% 
%     rightVectors = pos - rightCorr;
%     leftVectors = pos - leftCorr;
%     rightNormalsCorr = rightNormals(:, Idright');
%     leftNormalsCorr = leftNormals(:, Idleft');
% 
%     z1 = -(dot(rightVectors, rightNormalsCorr, 1) .* dot(leftVectors, leftNormalsCorr, 1)); % * 10^6;
    
    % N = size(TestTrack.cline, 2);
    N = size(cline_interp, 2);
    %xy_rep = repmat([x;y], 1, N);
%     [~,left_bound_ind] = min(norm(TestTrack.bl-xy_rep));
%     left_bound = TestTrack.bl(left_bound_ind);
%     [~,right_bound_ind] = min(norm(TestTrack.br-xy_rep));
%     right_bound = TestTrack.br(right_bound_ind);
    diff_norms = zeros(1,N);
    for i = 1:N
        diff_norms(i) = norm(cline_interp(:,i)-[x;y]);
    end
    [~,cline_ind] = min(diff_norms);
%     cpoint = TestTrack.cline(cline_ind+1);
%     next_cpoint = TestTrack.cline(cline_ind + 2);
%     next_next_cpoint = TestTrack.cline(cline_ind + 3);

     z2 = (N - cline_ind); % / N;  %;
     progress_bar = 100*(cline_ind / N)
     
%      z3 = norm(cline_interp(:,cline_ind+100) - [x;y] );
    
    % z2 = 2*norm([x;y] - cpoint)^2 + 5*norm([x;y] - next_cpoint)^2 + 5*norm([x;y] - next_next_cpoint)^2;
    % z2 = sqrt(z2)/10;% /10;

    % z = z1 + z2;
    %z = 10000*z1 + z2 + z3*10000; % was times
    % z = 10000000*z1 + z2 + z3*100000; % was times
    
    closest_center = cline_interp(:,cline_ind);
    z4 = norm([x;y] - closest_center)
    
    
    
    % z = z2; % + z4/2000000000000; %/20;
    
    z = z2 + 900*z4; % 500 okayish, 700 better
end

function [Lnormals, Rnormals] = calcNormals(leftBound, rightBound)
    n = size(leftBound, 2); % 246
    dLeft = leftBound(:,2:n) - leftBound(:,1:n-1);
    dRight = rightBound(:,2:n) - rightBound(:,1:n-1);
    Lnormals = [0, -1; 1, 0] * (dLeft ./ vecnorm(dLeft));
    Lnormals(:,n) = Lnormals(:,n - 1);
    Rnormals = [0, 1; -1, 0] * (dRight ./ vecnorm(dRight));
    Rnormals(:,n) = Rnormals(:,n - 1);
end

function z = hill(x,y, TestTrack, cline_interp)
    P = [x;y];
    % Can pass in interpolated TestTrack
    N = size(cline_interp, 2);
    diff_norms = zeros(1,N);
    for i = 1:N
        diff_norms(i) = norm(cline_interp(:,i)-P);
    end
    [~,cline_ind] = min(diff_norms);
    progress_bar = 100*(cline_ind / N)
    
    cpoint = cline_interp(:,cline_ind+1);    
    next_cpoint = cline_interp(:,cline_ind+5);
    next_next_cpoint = cline_interp(:,cline_ind+10);
    
    z = 1*norm(P - cpoint)^2 + 25*norm(P - next_cpoint)^2 + 50*norm(P - next_next_cpoint)^2;
    
    N = size(TestTrack.cline, 2);
    diff_norms_bl = zeros(1,N);
    diff_norms_br = zeros(1,N);
    for i = 1:N
        diff_norms_bl(i) = norm(TestTrack.bl(:,i)-P);
        diff_norms_br(i) = norm(TestTrack.br(:,i)-P);
    end
    [~,lb_ind] = min(diff_norms_bl);
    [~,rb_ind] = min(diff_norms_br);
    
    lb = TestTrack.bl(:,lb_ind);
    rb = TestTrack.br(:,rb_ind);
    
%     if norm(P - lb) > 10*norm(lb - cpoint)
%         disp("Outside barrier Left")
%         z = z + 10000;
%     elseif norm(P - rb) > 10*norm(rb - cpoint)
%         z = z + 10000;  
%         disp("Outside barrier Right")
%     end
    
    z;
    
    %z = inequalityConstraintTrack([x;y], left_bound, right_bound)
    %     xy_rep = repmat([x;y], 1, N);
%     [~,left_bound_ind] = min(norm(TestTrack.bl-xy_rep));
%     left_bound = TestTrack.bl(left_bound_ind);
%     [~,right_bound_ind] = min(norm(TestTrack.br-xy_rep));
%     right_bound = TestTrack.br(right_bound_ind);
%     [~,cline_ind] = min(norm(TestTrack.cline-xy_rep));
end

function dg = gradientHill(x,y,TestTrack,cline_interp,step_size)
    %f = @(x,y) hill(x,y,TestTrack,cline_interp) ;
    f = @(x,y) hill2(x,y,TestTrack,cline_interp) ;
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

% function dg = DinequalityConstraintTrack(pos, leftBound, rightBound)
%     dg = torGradient(inequalityConstraintTrack, pos, leftBound, rightBound);
% end

