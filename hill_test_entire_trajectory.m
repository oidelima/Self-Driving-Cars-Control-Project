clear
load TestTrack.mat
num_passed = 0;

% Interpolate Test Track?
num_interpl = 50; %10; % 5000; % 100 % 1000
cline_interp = interp_cline(TestTrack.cline, num_interpl);
% plot(TestTrack.cline(1,:), TestTrack.cline(2,:) )
% hold on
% plot(cline_interp(1,:), cline_interp(2,:) )
hillMap = makeHill(cline_interp);

%%
start_pos = [287,5,-176,0,2,0];
start_input = [0.01, 0.01];

Y(1,:) = start_pos;
U(1,:) = start_input;


i = 1;
while ~close_to_end(Y(i,:),TestTrack)
    % Curr pos
    curr_state = Y(i,:);
    curr_input = U(i,:);
    
    % Generate N inputs based on current position
    N = 3;
    U_temp = getControlInput(curr_state, curr_input,hillMap, N); % TestTrack, cline_interp, N);
    
    %Simulate forwards to get next state from U
    [Y_temp,T]=forwardIntegrateControlInput(U_temp,curr_state);
    
    
    Y(i+1:i+N,:) = Y_temp;
    U(i+1:i+N,:) = U_temp;
    i = i + N;
    
    if num_passed > 4
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

function U = getControlInput(curr_state, curr_input,hillMap,N) % TestTrack,cline_interp,N)
    % Calculate loss and derivative of loss with respect to input
    [L, dL_dU] = hilllossfun(curr_state, curr_input,hillMap);
    L;
    
    % Learning rate. Can tune
    eta = 1; % .3; % 3 best; %.5 % .3; % .1;% .3; % 10; % 1.8; % .5 gets closer to end % 300 still slowly gets better
    
    % Gradient Descent N steps, all the same
     U = ones(N,1) * (curr_input - eta * dL_dU);
    
    % Simulate N times in that direction
%     U(1,:) = (curr_input - eta * dL_dU);
%     for i = 1:N-1
%         U(i+1,:) = U(i,:) - eta * dL_dU;
%     end

end

function is_close = close_to_end(curr_state,TestTrack)
    x = curr_state(1)
    y = curr_state(3)
    v = curr_state(2)
    last_point = TestTrack.cline(:,end);
    dist = norm(last_point-[x;y]);
    is_close = ( dist < 10 );
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

function hillMap = makeHill(cline_interp)
    
    % ylim([-200,800])
    % xlim([200,1500])
    height = 819 - -200;
    width = 1500 - 200;
    hillMap = 5000000*ones(height,width);
    
    N = size(cline_interp,2); 
    
    for i = 1:N
        curr_x = cline_interp(1,i);
        curr_y = cline_interp(2,i);
        
        x = round(curr_x) - 199;
        y = round(curr_y) + 201;
        hillMap(y,x) = 5*(N-i);
       
        % road radius = 6
        for j = -5:5
            j_val = j;
            if y+j > 1000
                j_val = 1000-y;
            elseif y+j < 0
                j_val = -y+1;  
            end
            
            for k = -5:5
                k_val = k;
                if x+k > 1300
                    k_val = 1300-x;
                elseif x+k < 0
                    k_val = -x+1;  
                end
              
                if hillMap(y+j_val,x+k_val) > 5*(N-i)+4*(abs(k_val)+abs(j_val)) % N-i+k_val+j_val
                    hillMap(y+j_val,x+k_val) = 5*(N-i)+4*abs((k_val)+abs(j_val));
                end
            end
        end


    end
    

end


