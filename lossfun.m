function [L, dL_dU] = lossfun(curr_state, curr_input) %, goal_x, goal_y)
    
    % Step 1 forward integrate curr_input to get numerical gradient of
    % dx/dU and dy/dU with central difference
    step_size = .001;
    num_simulated = 5;
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
    dx_du1 = (cd1_dx_du1 + cd2_dx_du1) / 2
    
    % dx_du2
    cd1_dx_du2 = (x_u2p1 - x_0) / step_size;
    cd2_dx_du2 = (x_0 - x_u2m1) / step_size;
    dx_du2 = (cd1_dx_du2 + cd2_dx_du2) / 2
    
    % dy_du1
    cd1_dy_du1 = (y_u1p1 - y_0) / step_size;
    cd2_dy_du1 = (y_0 - y_u1m1) / step_size;
    dy_du1 = (cd1_dy_du1 + cd2_dy_du1) / 2
    
    % dy_du2
    cd1_dy_du2 = (y_u2p1 - y_0) / step_size;
    cd2_dy_du2 = (y_0 - y_u2m1) / step_size;
    dy_du2 = (cd1_dy_du2 + cd2_dy_du2) / 2
    
    
    
    
    
    
    
%     % Set up function forward integrate as only function of u1,u2
%     f_vec = @(u1,u2) forwardIntegrateControlInput([u1, u2], curr_state)
% 
% 
%     
%     
%     df_dx = brianNumGradient(f_x, curr_input(1), curr_input(2), step_size)
    




% temp_cell = num2cell(curr_state);
% [x_sym, u_sym, y_sym, v_sym, psi_sym, r_sym] = deal(temp_cell{:});
% 
% F_total = ((25952969516495739565766447169169*sin((6*atan((8*atan((243*atan2(v_sym - (29*r_sym)/20, u_sym))/(5*pi)))/5 - (3159*atan2(v_sym - (29*r_sym)/20, u_sym))/(25*pi)))/5)^2)/1208925819614629174706176 + 4*F_x_sym^2)^(1/2)
% F_max = 9.609879999999997e+03
% 
% if F_total <= F_max
% 
% % if F_total <= F_max
% %     L = (F_x_sym/70000 - y_goal_sym + y_sym + (r_sym*v_sym)/100 + (2735884546825257*sin((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*sin(delta_f_sym))/76965813944320000 - 4903/5000000)^2 + (x_goal_sym - x_sym - (u_sym*cos(psi_sym))/100 + (v_sym*sin(psi_sym))/100)^2
% %     dL_dU = [ 2*((2735884546825257*sin((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*cos(delta_f_sym))/76965813944320000 + (8207653640475771*cos((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*sin(delta_f_sym)*(1944/(25*pi*((729*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym))^2)/(10000*pi^2) + 1)) - 3159/(25*pi)))/(192414534860800000*(((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi))^2 + 1)))*(F_x_sym/70000 - y_goal_sym + y_sym + (r_sym*v_sym)/100 + (2735884546825257*sin((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*sin(delta_f_sym))/76965813944320000 - 4903/5000000), F_x_sym/2450000000 - y_goal_sym/35000 + y_sym/35000 + (r_sym*v_sym)/3500000 + (2735884546825257*sin((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*sin(delta_f_sym))/2693803488051200000000 - 4903/175000000000]
% % else
% %     L = (F_x_sym/70000 - y_goal_sym + y_sym + (r_sym*v_sym)/100 + (2735884546825257*sin((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*sin(delta_f_sym))/76965813944320000 - 4903/5000000)^2 + (x_goal_sym - x_sym - (u_sym*cos(psi_sym))/100 + (v_sym*sin(psi_sym))/100)^2
% %     dL_dU = [ 2*((2735884546825257*sin((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*cos(delta_f_sym))/76965813944320000 + (8207653640475771*cos((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*sin(delta_f_sym)*(1944/(25*pi*((729*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym))^2)/(10000*pi^2) + 1)) - 3159/(25*pi)))/(192414534860800000*(((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi))^2 + 1)))*(F_x_sym/70000 - y_goal_sym + y_sym + (r_sym*v_sym)/100 + (2735884546825257*sin((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*sin(delta_f_sym))/76965813944320000 - 4903/5000000), F_x_sym/2450000000 - y_goal_sym/35000 + y_sym/35000 + (r_sym*v_sym)/3500000 + (2735884546825257*sin((6*atan((8*atan((27*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(100*pi)))/5 - (351*(180*delta_f_sym - 180*atan2((27*r_sym)/20 + v_sym, u_sym)))/(500*pi)))/5)*sin(delta_f_sym))/2693803488051200000000 - 4903/175000000000]
% % end

end