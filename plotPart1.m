global n input_cost precompFlag
setGlobaln(10)
niter = 2;
subdivision_num = 5;
gradient_precision = [0.05, 0.05];

xlim1 = 200;
xlim2 = 1500;
ylim1 = -200;
ylim2 = 825;

input_cost = 0;
precompFlag = 0;

z0 = repmat([287, 5, -176, 0, 2, 0], n, 1);
u0 = repmat([0 0], n-1, 1);

x0 = encodeColocationVector(z0(:,1), z0(:,2), z0(:,3), z0(:,4), z0(:,5), z0(:,6), u0(:,1), u0(:,2));

%[z_out, breakpts] = generateMultiplePart1(n, niter, subdivision_num, gradient_precision);
%[x, u, y, v, psi, r, delta, Fx] = decodeColocationVector(z_out);
tic
num_pts = 11800;
gains = [5, 10000;
    0,0;
    0,0];
speeds = [15,70];
lookahead = 2;
cline_strength = 0.1;
[Yplot, U] = PDControlPart1(num_pts, z0(1,:), gains, speeds, lookahead, cline_strength);
toc

x = Yplot(:,1);
y = Yplot(:,3);
x_g = Yplot(:,7);
y_g = Yplot(:,9);
delta = U(:,1);
Fx = U(:,2);

ROB535_ControlsProject_part1_input = [delta, Fx];
save('ROB535_ControlsProject_part1_Team1.mat','ROB535_ControlsProject_part1_input')
[sol1, T] = forwardIntegrateControlInput(ROB535_ControlsProject_part1_input);

%% Plotting

height_vis = 1:num_pts;
hold on

plot3(x_g,y_g, height_vis, 'o-g')
plot(x,y, 'o-b')
%plot(breakpts(:,1),breakpts(:,3), 'ob')
plot(sol1(:,1), sol1(:,3), '.r')

load TestTrack.mat
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), '--k')
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
xlim([xlim1,xlim2])
ylim([ylim1,ylim2])

figure(2)
hold on
plot(Yplot(:, 5), '.-b');
plot(Yplot(:, 5+6), '.-r');
plot(Yplot(:, 5+12), '.-k');
plot(Yplot(:, 2), 'o-b');
plot(Yplot(:, 2+6), 'o-r');
plot(Yplot(:, 2+12), 'o-k');
plot(Yplot(:, 1+18), '--m');