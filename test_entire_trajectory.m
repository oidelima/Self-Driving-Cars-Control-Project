clear
load TestTrack.mat

% Interpolate Test Track?
num_interpl = 5;
cline_interp = interp_cline(TestTrack.cline, num_interpl);
plot(TestTrack.cline(1,:), TestTrack.cline(2,:) )
hold on
plot(cline_interp(1,:), cline_interp(2,:) )

start_pos =[287,5,-176,0,2,0]';

Y(1,:) = start_pos;

% Window size when segmenting track
N = 100;




%for i = 1:N:size(cline_interp,2)-N+1
for i = 1
    % Start pos
    curr_state = Y(i,:);
    
    % Estimate the test_track
    %track_section = TestTrack.cline(:,i:i+5);
    track_section = cline_interp(:,i:i+N-1);
    
    %Generate N inputs based on current position
    U = getControlInput(curr_state, track_section);
    
    %Simulate forwards to get next state from U
    [Y_temp,T]=forwardIntegrateControlInput(U,curr_state);
    
    
    Y(i+1:i+N,:) = Y_temp;
    
%     subplot(2,1,1)
%     plot(track_section(1,:),track_section(2,:))
%     
%     subplot(2,1,2)
%     plot( Y(1,:), Y(3,:) )
end

function U = getControlInput(curr_state, track_section)
    U = zeros(size(track_section,2));
end


function cline_interp = interp_cline(cline, num_interpl)
    num_points = size(cline,2);
    
    for i = 1:num_points
        x = cline(1,i);
        y = cline(2,i);
    end
    
    cline_interp = cline;
end
