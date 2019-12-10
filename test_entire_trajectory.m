
load TestTrack.mat
load ROB535_ControlsProject_part1_Team1.mat

plot(TestTrack.cline(1,:), TestTrack.cline(2,:),'r' )
hold on
plot(TestTrack.br(1,:), TestTrack.br(2,:),'r' )
hold on
plot(TestTrack.bl(1,:), TestTrack.bl(2,:),'r' )
hold on
%Xobs = generateRandomObstacles(9 + randi(16),TestTrack);

[Y,U,Xobs, t_total,t_update] = forwardIntegrate();
info = getTrajectoryInfo(Y,U,Xobs,TestTrack)
plot(Y(:,1), Y(:,3),'k' )
%[Uref, Yref,T] = create_reference();  % when there are no obstacles 
%[Yref,U,t_total,t_update] = forwardIntegrate(); %when there are obstacles
%plot(Yref(:,1), Yref(:,3),'k')

for i = 1 : length(Xobs)
Obs= cell2mat(Xobs(i));

x = [Obs(1,1), Obs(2,1), Obs(3,1), Obs(4,1), Obs(1,1)];
y = [Obs(1,2), Obs(2,2), Obs(3,2), Obs(4,2), Obs(1,2)];
plot(x, y, 'b-', 'LineWidth', 0.5);
hold on
end



function [Uref, Yref,T] = create_reference()
    load ROB535_ControlsProject_part1_Team1.mat

    x0 =[287,5,-176,0,2,0];
    Uref = ROB535_ControlsProject_part1_input;
   [Yref,T]=forwardIntegrateControlInput(ROB535_ControlsProject_part1_input,x0);
end