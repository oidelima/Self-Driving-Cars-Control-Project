clear
global n input_cost
setGlobaln(10)
niter = 50;
subdivision_num = 6;
gradient_precision = [0.05, 0.5];

input_cost = 0;

z0 = repmat([287, 5, -176, 0, 2, 0], n, 1);
u0 = repmat([0 0], n-1, 1);

x0 = encodeColocationVector(z0(:,1), z0(:,2), z0(:,3), z0(:,4), z0(:,5), z0(:,6), u0(:,1), u0(:,2));

tic
[z_out, breakpts] = generateMultiplePart1(n, niter, subdivision_num, gradient_precision);
toc

% Testing inequality constraint fitter
xlim1 = 200;
xlim2 = 1500;
ylim1 = -200;
ylim2 = 825;
hold on

[x, u, y, v, psi, r, delta, Fx] = decodeColocationVector(z_out);
plot(x,y, '.-b')
plot(breakpts(:,1),breakpts(:,3), 'ob')

ROB535ControlsProjectpart1input = [delta, Fx];
[sol1, ~] = forwardIntegrateControlInput(ROB535ControlsProjectpart1input);
plot(sol1(:,1), sol1(:,3), '.r')

load TestTrack.mat
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), '--k')
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
xlim([xlim1,xlim2])
ylim([ylim1,ylim2])
%view(0,90)
%{
for i = -2:10
    intermediate_lines = (i*TestTrack.bl + (8-i)*TestTrack.br)/8;
    plot(intermediate_lines(1,:), intermediate_lines(2,:), '.b')
end
%}