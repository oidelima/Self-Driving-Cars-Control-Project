clear
load TestTrack.mat

% Interpolate Test Track?
% num_interpl = 5;
% cline_interp = interp_cline(TestTrack.cline, num_interpl);
% plot(TestTrack.cline(1,:), TestTrack.cline(2,:) )
% hold on
% plot(cline_interp(1,:), cline_interp(2,:) )

start_pos = [287,5,-176,0,2,0];
start_input = [0.01, 0.01];

Y(1,:) = start_pos;
U(1,:) = start_input;

% % Window size when segmenting track
% M = 100;


i = 1;
%for i = 1:N:size(cline_interp,2)-N+1
while ~close_to_end(Y(i,:),TestTrack)
    % Start pos
    curr_state = Y(i,:);
    curr_input = U(i,:);
    
%     % Estimate the test_track
%     %track_section = TestTrack.cline(:,i:i+5);
%     track_section = cline_interp(:,i:i+N-1);
    
    % Generate N inputs based on current position
    N = 3;
    U_temp = getControlInput(curr_state, curr_input,TestTrack, N);
    
    %Simulate forwards to get next state from U
    [Y_temp,T]=forwardIntegrateControlInput(U_temp,curr_state);
    
    
    Y(i+1:i+N,:) = Y_temp;
    U(i+1:i+N,:) = U_temp;
    i = i + N;
end

function U = getControlInput(curr_state, curr_input,TestTrack,N)
    % Calculate loss and derivative of loss with respect to input
    [L, dL_dU] = lossfun(curr_state, curr_input,TestTrack);
    L;
    
    % Learning rate. Can tune
    eta = .5; % .5 gets closer to end % 300 still slowly gets better
    
    % Gradient Descent N steps
    U = ones(N,1) * (curr_input - eta * dL_dU);
%     U(1,:) = (curr_input - eta * dL_dU);
%     for i = 1:N-1
%         U(i+1,:) = U(i,:) - eta * dL_dU;
%     end
end

function is_close = close_to_end(curr_state,TestTrack)
    x = curr_state(1);
    y = curr_state(3);
    last_point = TestTrack.cline(:,end);
    dist = norm(last_point-[x;y])
    is_close = ( dist < 20 );
end


% function cline_interp = interp_cline(cline, num_interpl)
%     num_points = size(cline,2);
%     
%     for i = 1:num_points
%         x = cline(1,i);
%         y = cline(2,i);
%     end
%     
%     cline_interp = cline;
% end


