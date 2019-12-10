xlim1 = 200;
xlim2 = 1500;
ylim1 = -200;
ylim2 = 825;

plot(Y(:,1), Y(:,3), '.r')
hold on
load TestTrack.mat

num_cells = numel(Xobs);

for i=1:num_cells
    obs = Xobs{i};
    obs = [obs;obs(1,:)];
    plot(obs(:,1),obs(:,2))
end
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), '--k')
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
xlim([xlim1,xlim2])
ylim([ylim1,ylim2])