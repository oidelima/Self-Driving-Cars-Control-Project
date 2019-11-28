load TestTrack.mat
global Jdata dJdata

rng(0,'twister');
xlim1 = 200;
xlim2 = 1500;
ylim1 = -200;
ylim2 = 825;

hold on
rx = reshape(Jdata{1}, 1,[]);
ry = reshape(Jdata{2}, 1,[]);
rz = reshape(Jdata{3}, 1,[]);
passmap = inequalityConstraintTrack([rx;ry], TestTrack.bl, TestTrack.br) <= 0;
plot3(passmap .* rx, passmap .* ry, passmap .* rz, '.b')

plot(TestTrack.cline(1,:), TestTrack.cline(2,:), '--k')
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
xlim([xlim1,xlim2])
ylim([ylim1,ylim2])
view(0,90)
%{
for i = -2:10
    intermediate_lines = (i*TestTrack.bl + (8-i)*TestTrack.br)/8;
    plot(intermediate_lines(1,:), intermediate_lines(2,:), '.b')
end
%}