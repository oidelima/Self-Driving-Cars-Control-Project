generatePart1

load TestTrack.mat

sol1 = forwardIntegrateControlInput(ROB535ControlsProjectpart1input);

plot(TestTrack.cline(1,:), TestTrack.cline(2,:), '--k')
hold on
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'ok')
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'ok')
plot(sol1(:,1), sol1(:,3), 'b')
%xlim([200,1500])
%ylim([-200, 825])