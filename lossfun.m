function [L, dL_dU] = lossfun(curr_state, curr_input, TestTrack) %, goal_x, goal_y)
    x = curr_state(1);
    y = curr_state(3);

    L = hill(x,y, TestTrack);
    
    hill_step = .05;
    dL_dxy = gradientHill(x,y,TestTrack,hill_step);

    % Step 1 forward integrate curr_input to get numerical gradient of
    % dx/dU and dy/dU with central difference
    step_size = .001;
    num_simulated = 5; % >= 3 or the program crashes in ode45
    u1 = curr_input(1);
    u2 = curr_input(2);
    
    % Simulate the resulting Y from the inputs
    [output_0,~] = forwardIntegrateControlInput([u1*ones(num_simulated,1), u2*ones(num_simulated,1)], curr_state);
    
    % Shift u1 by step size and forward integrate
    [output_u1p1,~] = forwardIntegrateControlInput([(u1+step_size)*ones(num_simulated,1), u2*ones(num_simulated,1)], curr_state);
    [output_u1m1,~] = forwardIntegrateControlInput([(u1-step_size)*ones(num_simulated,1), u2*ones(num_simulated,1)], curr_state);
    
    % Shift u2 by step size and forward integrate
    [output_u2p1,~] = forwardIntegrateControlInput([u1*ones(num_simulated,1), (u2+step_size)*ones(num_simulated,1)], curr_state);
    [output_u2m1,~] = forwardIntegrateControlInput([u1*ones(num_simulated,1), (u2-step_size)*ones(num_simulated,1)], curr_state);
    
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
    cd1_dx_du1 = (x_u1p1 - x_0) / step_size;
    cd2_dx_du1 = (x_0 - x_u1m1) / step_size;
    dx_du1 = (cd1_dx_du1 + cd2_dx_du1) / 2;
    
    % dx_du2
    cd1_dx_du2 = (x_u2p1 - x_0) / step_size;
    cd2_dx_du2 = (x_0 - x_u2m1) / step_size;
    dx_du2 = (cd1_dx_du2 + cd2_dx_du2) / 2;
    
    % dy_du1
    cd1_dy_du1 = (y_u1p1 - y_0) / step_size;
    cd2_dy_du1 = (y_0 - y_u1m1) / step_size;
    dy_du1 = (cd1_dy_du1 + cd2_dy_du1) / 2;
    
    % dy_du2
    cd1_dy_du2 = (y_u2p1 - y_0) / step_size;
    cd2_dy_du2 = (y_0 - y_u2m1) / step_size;
    dy_du2 = (cd1_dy_du2 + cd2_dy_du2) / 2;
    
    
    % dL_dU
    dL_dU = dL_dxy * [dx_du1 dx_du2; dy_du1 dy_du2];
    
  

end

% function g = inequalityConstraintTrack(pos, leftBound, rightBound)
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
%     g = -(dot(rightVectors, rightNormalsCorr, 1) .* dot(leftVectors, leftNormalsCorr, 1));
% end

function [Lnormals, Rnormals] = calcNormals(leftBound, rightBound)
    n = size(leftBound, 2); % 246
    dLeft = leftBound(:,2:n) - leftBound(:,1:n-1);
    dRight = rightBound(:,2:n) - rightBound(:,1:n-1);
    Lnormals = [0, -1; 1, 0] * (dLeft ./ vecnorm(dLeft));
    Lnormals(:,n) = Lnormals(:,n - 1);
    Rnormals = [0, 1; -1, 0] * (dRight ./ vecnorm(dRight));
    Rnormals(:,n) = Rnormals(:,n - 1);
end

function z = hill(x,y, TestTrack)
    % Can pass in interpolated TestTrack
    N = size(TestTrack, 2);
    xy_rep = repmat([x;y], 1, N);
    [~,left_bound_ind] = min(norm(TestTrack.bl-xy_rep));
    left_bound = TestTrack.bl(left_bound_ind);
    [~,right_bound_ind] = min(norm(TestTrack.br-xy_rep));
    right_bound = TestTrack.br(right_bound_ind);
    [~,cline_ind] = min(norm(TestTrack.cline-xy_rep));
    cpoint = TestTrack.cline(cline_ind);
    
    next_cpoint = TestTrack.cline(cline_ind + 1);
    next_next_cpoint = TestTrack.cline(cline_ind + 2);
    
    z = 2*norm([x;y] - cpoint)^2 + 5*norm([x;y] - next_cpoint)^2 + norm([x;y] - next_next_cpoint)^2;
    
    if norm([x;y] - cpoint) > norm(cpoint - left_bound)
        z = z + 10000;
    elseif norm([x;y] - cpoint) > norm(cpoint - right_bound)
       z = z + 10000;     
    end
    
    %z = inequalityConstraintTrack([x;y], left_bound, right_bound)
end

function dg = gradientHill(x,y,TestTrack, step_size)
    f = @(x,y) hill(x,y,TestTrack) ;
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

