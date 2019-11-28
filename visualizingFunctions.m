clear
load TestTrack.mat

rng(0,'twister');
xlim1 = 200;
xlim2 = 1500;
ylim1 = -200;
ylim2 = 825;

new_cline = subdivideTrack(TestTrack.cline);
new_cline = subdivideTrack(new_cline);
new_cline = subdivideTrack(new_cline);

rx = (xlim2-xlim1).*rand(100000,1) + xlim1;
ry = (ylim2-ylim1).*rand(100000,1) + ylim1;
j = costFunctionTrack([rx';ry'], new_cline);

hold on
plot3(rx, ry, j, '.b')

plot(TestTrack.cline(1,:), TestTrack.cline(2,:), 'o--k')
plot(new_cline(1,:), new_cline(2,:), 'o--r')
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
xlim([xlim1,xlim2])
ylim([ylim1,ylim2])
view(0,90)