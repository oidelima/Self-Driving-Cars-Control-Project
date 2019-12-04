load TestTrack.mat
%plot(Y(:,1),Y(:,3))
% hold on
% plot(TestTrack.cline(1,:),TestTrack.cline(2,:))
% plot(TestTrack.br(1,:),TestTrack.br(2,:))
% plot(TestTrack.bl(1,:),TestTrack.bl(2,:))
% ylim([-200,800])
% xlim([200,1500])

% Interpolate Test Track?
num_interpl = 10 % 5000; % 100 % 1000
cline_interp = interp_cline(TestTrack.cline, num_interpl);

hillMap = makeHill(cline_interp);
surf(hillMap)

hill(283,300, hillMap)

function z = hill(x,y,hillMap)
    if y < -200 || y > 800 || x < 200 || x > 1500
        z = 1000;
    else
        x = round(x);
        y = round(y);
        z = hillMap(y+201,x-199);
    end
end

function hillMap = makeHill(cline_interp)
    
    % ylim([-200,800])
    % xlim([200,1500])
    height = 819 - -200;
    width = 1500 - 200;
    hillMap = 50000*ones(height,width);
    
    N = size(cline_interp,2); 
    
    for i = 1:N-1
        curr_x = cline_interp(1,i);
        curr_y = cline_interp(2,i);
        
        x = round(curr_x) - 199;
        y = round(curr_y) + 201;
        hillMap(y,x) = 5*(N-i);

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
              
                if hillMap(y+j_val,x+k_val) > 5*(N-i)+3*abs(k_val+j_val) % N-i+k_val+j_val
                    hillMap(y+j_val,x+k_val) = 5*(N-i)+3*abs(k_val+j_val);
                end
            end
        end


    end
    
   
    % road radius = 6
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
