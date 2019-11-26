clear
load TestTrack.mat
num_passed = 0

% Interpolate Test Track?
num_interpl = 100; % 1000
cline_interp = interp_cline(TestTrack.cline, num_interpl);
% plot(TestTrack.cline(1,:), TestTrack.cline(2,:) )
% hold on
% plot(cline_interp(1,:), cline_interp(2,:) )

%%
start_pos = [287,5,-176,0,2,0];
start_input = [0.01, 0.01];

Y(1,:) = start_pos;
U(1,:) = start_input;

% % Window size when segmenting track
% M = 100;


i = 1;
%for i = 1:N:size(cline_interp,2)-N+1
while ~close_to_end(Y(i,:),TestTrack)
    % Curr pos
    curr_state = Y(i,:);
    curr_input = U(i,:);
    
    % Generate N inputs based on current position
    N = 3;
    U_temp = getControlInput(curr_state, curr_input,TestTrack, cline_interp, N);
    
    %Simulate forwards to get next state from U
    [Y_temp,T]=forwardIntegrateControlInput(U_temp,curr_state);
    
    
    Y(i+1:i+N,:) = Y_temp;
    U(i+1:i+N,:) = U_temp;
    i = i + N;
    
    if num_passed > 2
        close all;
        plot(Y(:,1),Y(:,3))
        hold on
        plot(TestTrack.cline(1,:),TestTrack.cline(2,:))
        plot(TestTrack.br(1,:),TestTrack.br(2,:))
        plot(TestTrack.bl(1,:),TestTrack.bl(2,:))
        ylim([-200,800])
        xlim([200,1500])
        num_passed = 0;
    end
    num_passed = num_passed + 1;
    
end

function U = getControlInput(curr_state, curr_input,TestTrack,cline_interp,N)
    % Calculate loss and derivative of loss with respect to input
    [L, dL_dU] = lossfun(curr_state, curr_input,TestTrack,cline_interp);
    L;
    
    % Learning rate. Can tune
    eta = .3; % 10; % 1.8; % .5 gets closer to end % 300 still slowly gets better
    
    % Gradient Descent N steps, all the same
    U = ones(N,1) * (curr_input - eta * dL_dU);
    
    % Simulate N times in that direction
%     U(1,:) = (curr_input - eta * dL_dU);
%     for i = 1:N-1
%         U(i+1,:) = U(i,:) - eta * dL_dU;
%     end

end

function is_close = close_to_end(curr_state,TestTrack)
    x = curr_state(1);
    y = curr_state(3);
    last_point = TestTrack.cline(:,end);
    dist = norm(last_point-[x;y]);
    is_close = ( dist < 20 );
end


function cline_interp = interp_cline(cline, num_interpl)
    N = size(cline,2);
    step_size = 1/num_interpl;
    
    x_points = cline(1,:);
    y_points = cline(2,:);
    new_ind = 1:step_size:N;
    
    x_interp = interp1(x_points, new_ind);
    y_interp = interp1(y_points, new_ind);
    
    cline_interp = [x_interp;y_interp];
end


