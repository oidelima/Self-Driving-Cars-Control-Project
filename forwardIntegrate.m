
function [Y,U,Xobs, t_total,t_update] = forwardIntegrate()
% [Y,U,t_total,t_update] = forwardIntegrate
% 
% This script returns the vehicle trajectory with control input being
% generated via the control input generation function:
%         ROB535_ControlsProject_part2_Team<your team number>
% Obstacles are randomly generated along the test track. Notice that the
% vehicle can only sense (observe) the obstacles within 150m, therefore
% the control input generation function is called repeatedly. In this
% script, we assume the control input generation function is called every
% 'dt' second (see line 32). 
% 
% OUTPUTS:
%   Y         an N-by-6 vector where each column is the trajectory of the
%             state of the vehicle 
%   
%   U         an N-by-2 vector of inputs, where the first column is the
%             steering input in radians, and the second column is the
%             longitudinal force in Newtons
%   
%   t_total   a scalar that records the total computational time  
%   
%   t_update  a M-by-1 vector of time that records the time consumption
%             when the control input generation function is called
%             
% Written by: Jinsun Liu
% Created: 31 Oct 2019


    load('TestTrack.mat') % load test track

    dt = 0.5; 
    TOTAL_TIME = 20*60; % second
    
    % initialization
    t_total = 0;
    t_update = zeros(TOTAL_TIME/dt+1,1);
    Y = zeros(TOTAL_TIME/0.01+1,6);
    U = zeros(TOTAL_TIME/0.01,2);
    Y(1,:) = [287,5,-176,0,2,0];

    % generate obstacles along the track
    %Xobs = generateRandomObstacles(9 + randi(16),TestTrack);

    %Xobs = {[252.123205149458,-76.8012453452512;255.813628930957,-75.5052510460710;254.990072825209,-73.1601241886638;251.299649043709,-74.4561184878440],[245.122040752708,46.4984881232945;249.672147879113,44.4441378869165;250.317841299765,45.8742612333431;245.767734173360,47.9286114697211],[340.438067835437,116.735711621455;343.816302786632,111.061863576539;345.951504393059,112.333172331651;342.573269441864,118.007020376567],[370.036525593549,141.327877425725;372.870163281737,136.850586394358;374.089332727349,137.622187897403;371.255695039162,142.099478928771],[467.019261807555,221.169254250867;472.945906912304,218.805280292571;473.377496050605,219.887303909961;467.450850945855,222.251277868258],[564.721992952314,320.028876782073;569.097481134199,315.990077223612;571.506789406232,318.600233959099;567.131301224348,322.639033517560],[585.025611776273,397.620690281571;590.189747871196,396.246106653173;590.879535903692,398.837552675711;585.715399808769,400.212136304110],[689.330103542752,474.954547961345;689.805165033862,470.996013886167;693.566999005468,471.447469506332;693.091937514359,475.406003581510],[779.332504536669,451.322729944863;777.072187473337,446.127470068381;779.905634510850,444.894713894637;782.165951574182,450.089973771119],[810.331530883775,431.058218740067;808.917313363906,428.039098600750;811.406167310208,426.873268593054;812.820384830077,429.892388732371],[849.690711860756,431.735329492723;852.761539924372,427.167807156680;855.621463127468,429.090585368961;852.550635063852,433.658107705005],[902.529835417303,476.789110872205;906.875757464676,473.188389515618;909.196184178093,475.989048776493;904.850262130720,479.589770133079],[995.967135302213,515.136193688979;993.163978103756,509.923041228732;996.657670890647,508.044452378809;999.460828089104,513.257604839056],[1099.09978561681,483.488157378343;1100.50782657697,477.077970113517;1101.48607197622,477.292848327573;1100.07803101606,483.703035592399],[1216.73104127047,509.098003375503;1218.24100236674,502.335310954287;1221.75089271119,503.118992592038;1220.24093161492,509.881685013255],[1259.07024033186,546.921211091997;1263.81086891666,543.120470320181;1265.58591061006,545.334463301065;1260.84528202526,549.135204072881],[1314.11954427958,613.775836904778;1318.22209669704,610.379801029937;1320.75350806736,613.437851127060;1316.65095564990,616.833887001901],[1332.03888912267,635.245158413063;1337.17678518034,630.924304346702;1338.83930870425,632.901198943955;1333.70141264659,637.222053010315],[1347.44559156352,662.457002815341;1351.57576444419,659.103211218799;1353.11296377780,660.996262357222;1348.98279089712,664.350053953764],[1373.17014625668,685.466918171726;1375.51190701104,683.484908632757;1377.71105805999,686.083223887091;1375.36929730563,688.065233426060]}
    Xobs =  {[231.939167454530,-26.6595863430617;236.050553840232,-25.7761764436145;235.449136939807,-22.9771844263346;231.337750554105,-23.8605943257817],[260.039130672869,66.3258392803640;263.437834492277,63.8888189003786;265.023495168351,66.1002041806650;261.624791348943,68.5372245606504],[345.523680321902,124.865834684383;347.730500774146,121.138374511580;348.725438110456,121.727421197710;346.518617658213,125.454881370513],[422.505853545132,167.243674503119;425.070620828948,163.105794812999;428.295065800020,165.104391096505;425.730298516203,169.242270786625],[450.734706724693,190.608009960173;454.963671653253,186.948065128004;457.048592049210,189.357132924831;452.819627120649,193.017077757000],[494.376407316321,274.245774319094;496.641721030210,272.128225379456;498.836170528096,274.475805852779;496.570856814206,276.593354792416],[580.314273314352,379.489530012441;586.634341222271,378.702960253671;586.781336269884,379.884061704228;580.461268361966,380.670631462998],[694.147066107716,480.500517398323;694.495198562488,475.275973106677;698.068595927405,475.514082983267;697.720463472633,480.738627274913],[751.177413623902,457.879739278949;748.401634176218,452.094824861781;751.830546969276,450.449527287347;754.606326416961,456.234441704515],[797.086445971437,437.355442505388;795.299069324644,433.344120487912;798.910017134519,431.735143756540;800.697393781312,435.746465774016],[889.065702975230,461.924708777929;892.384511483340,458.478959806772;894.960710689023,460.960252528318;891.641902180913,464.406001499475],[981.274446087659,521.747076682878;980.444246615084,518.639089036654;983.687220489813,517.772832311456;984.517419962388,520.880819957680],[1042.83196120780,475.284906000886;1042.51356462974,471.590169516857;1045.04431726025,471.372080085862;1045.36271383832,475.066816569891],[1091.33000428216,481.845275713152;1092.62287849549,476.343921135634;1094.11503392902,476.694592821145;1092.82215971569,482.195947398663],[1247.05869493673,532.486185612121;1250.48238818945,529.846226222050;1251.82319688117,531.585085421291;1248.39950362845,534.225044811362],[1266.33724681939,555.764318929183;1271.43342615292,551.579935946394;1272.32297593296,552.663322726205;1267.22679659944,556.847705708993],[1285.21253804976,584.260559898008;1287.88434617488,582.042331659084;1288.58231712846,582.883022635854;1285.91050900334,585.101250874777],[1307.85147852355,617.111299198392;1313.21095549315,612.680773616226;1315.64322442057,615.623018065253;1310.28374745097,620.053543647419]}
    %Xobs = {[243.416218619381,29.3837751452655;248.254617066221,27.6635928130948;248.683436466143,28.8697435098874;243.845038019303,30.5899258420580],[261.538693245829,61.7360853326646;264.391857648953,59.8204161937293;265.463291897583,61.4161917879314;262.610127494458,63.3318609268667],[385.650781143623,151.954297512734;388.961945915595,146.590813126927;390.977047374112,147.834842657255;387.665882602140,153.198327043061],[482.848241337904,258.728661512938;486.599921488372,256.522261191331;487.261482123098,257.647154021217;483.509801972631,259.853554342824],[572.856253343147,330.204575347809;576.547074770188,328.251139154631;578.411760628593,331.774275740784;574.720939201551,333.727711933962],[580.715763080056,379.585336423216;586.654482574355,378.840438610859;586.902930436320,380.821196423130;580.964210942021,381.566094235487],[666.174496746423,469.461886354033;667.360382392161,466.661267921198;669.102664437231,467.399014714496;667.916778791494,470.199633147331],[736.354219580690,464.822907099718;735.057280187592,461.904304922716;736.141214581527,461.422636942713;737.438153974625,464.341239119715],[777.476220116344,452.085278266311;775.173984516203,446.901892376129;776.905597290735,446.132784951050;779.207832890876,451.316170841232],[827.526923335289,424.608803273668;826.384086717354,419.316928955736;827.497472773563,419.076481379227;828.640309391498,424.368355697159],[969.572257013870,529.102603714063;969.566401982135,523.998197208899;972.082684754129,523.995310895775;972.088539785865,529.099717400939],[1001.90866770884,518.235226099350;998.670498835932,513.280742960908;1001.41967001753,511.483929792350;1004.65783889044,516.438412930792],[1041.84637548824,477.146624144806;1040.97051803873,471.383788275588;1043.83088078955,470.949059596778;1044.70673823906,476.711895465996],[1068.61040305897,476.792789309868;1070.00406639038,470.827929435980;1071.24044269090,471.116803332482;1069.84677935949,477.081663206371],[1164.83390868346,497.968535595719;1165.77758887897,493.925293466055;1166.98765454406,494.207719046019;1166.04397434855,498.250961175684],[1185.24161252630,502.242647209694;1186.39958222020,496.589906043147;1189.66352369188,497.258527725178;1188.50555399798,502.911268891725],[1312.73616898746,620.298555076014;1316.75917022514,616.961844425789;1319.13108326341,619.821609652350;1315.10808202573,623.158320302575]}
    
    iteration = 1; % a counter that counts how many times the control input 
                   % generation function is called.

    TIMER = tic; % start the timer
    
    % you only have TOTAL_TIME seconds to sense the obstacles, update
    % control inputs, and simulate forward vehicle dynamcis.
    while t_total < TOTAL_TIME 
        curr_pos = Y( (iteration-1)*dt/0.01+1 , [1,3] ); % record current vehicle position
        Xobs_seen = senseObstacles(curr_pos, Xobs); % sense the obstacles within 150m
        curr_state = Y( (iteration-1)*dt/0.01+1 , : ); % record current vehicle states - 1 : 6 vector

        % compute control inputs, and record the time consumption
        t_temp = toc(TIMER);
        %%%%%%%%%%%%%%%% THIS IS WHERE YOUR FUNCTION IS CALLED (replace in your team number). %%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Utemp, FLAG_terminate] = ROB535_ControlsProject_part2_Team1(TestTrack,Xobs_seen,curr_state); %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %FLAG_terminate = randi(2)-1                                 % GSIs: This line is just for us to debug. Feel free to play with it if you want
%         Utemp = rand(dt/0.01+1 + FLAG_terminate* (randi(10)-5),2);  % GSIs: This line is just for us to debug. Feel free to play with it if you want        
        t_update(iteration) = toc(TIMER)-t_temp;
        
        
        
        % Utemp must contain control inputs for at least dt second,
        % otherwise stop the whole computation.
        if size(Utemp,1)<dt/0.01+1 && FLAG_terminate == 0
            fprintf('When FLAG_terminate = 0, Utemp cannot contain control inputs for less than %f second. \n',dt);
            fprintf('Solving process is terminated.\n');
            t_total = toc(TIMER);
            break
        end

        
        
        if FLAG_terminate == 0
            % if FLAG_terminate == 0, simulate forward vehicle dynamics for
            % dt second.
            U( (iteration-1)*dt/0.01+1:iteration*dt/0.01 , : ) = Utemp(1:dt/0.01,:);
            Ytemp = forwardIntegrateControlInput( Utemp(1:dt/0.01+1,:) , curr_state );
            Y( (iteration-1)*dt/0.01+2:iteration*dt/0.01+1 , : ) = Ytemp(2:end,:);
        
            % update the counter
            iteration = iteration + 1;
        else
            % if FLAG_terminate == 1, simulate forward vehicle dynamics for
            % no more than dt second, and stop the solving process.
            simulate_length = min(dt/0.01+1, size(Utemp,1));
            U((iteration-1)*dt/0.01+1:(iteration-1)*dt/0.01+simulate_length-1, :) = Utemp(1:simulate_length-1,:);
            Ytemp = forwardIntegrateControlInput( Utemp(1:simulate_length,:) , curr_state );
            Y((iteration-1)*dt/0.01+2:(iteration-1)*dt/0.01+simulate_length, : ) = Ytemp(2:end,:);
        end
        
        
        % update t_total
        t_total = toc(TIMER);

        % stop the computation if FLAG_terminate == 1
        if FLAG_terminate == 1
            break
        end
    end

    % if reach the finish line before TOTAL_TIME, ignore any parts of the
    % trajectory after crossing the finish line.
    idx_temp = find(sum(abs(Y),2)==0,1);
    Y(idx_temp:end,:) = [];
    U(idx_temp:end,:) = [];
    t_update(iteration:end) = [];

end