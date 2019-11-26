clear
generatePart1

load TestTrack.mat

%sol1 = forwardIntegrateControlInput(ROB535ControlsProjectpart1input);

% Testing inequality constraint fitter
rng(0,'twister');
%xlim1 = 600;
%xlim2 = 1100;
xlim1 = 200;
xlim2 = 1500;
%rx = (xlim2-xlim1).*rand(100000,1) + xlim1;
ylim1 = -200;
ylim2 = 825;
%ylim1 = 300;
%ylim2 = 600;
%ry = (ylim2-ylim1).*rand(100000,1) + ylim1;
[rx, ry] = meshgrid(xlim1:1:xlim2,ylim1:1:ylim2);
rx2 = reshape(rx,[1,size(rx,1)*size(rx,2)]);
ry2 = reshape(ry,[1,size(ry,1)*size(ry,2)]);
g2 = inequalityConstraintTrack([rx2;ry2], TestTrack.bl, TestTrack.br);
g3 = inequalityConstraintTrack(TestTrack.cline, TestTrack.bl, TestTrack.br);
g4 = inequalityConstraintTrack(TestTrack.cline + [1;1], TestTrack.bl, TestTrack.br);
dg2 = torGradient(@inequalityConstraintTrack, [rx2;ry2], 0.05, TestTrack.bl, TestTrack.br);
%passmap = g2 < 0;
%failmap = g2 > 0;
%passpoints = [passmap .* rx'; passmap .* ry'; passmap .* g2];
%passpoints(:,~any(passpoints))=[];
%failpoints = [failmap .* rx'; failmap .* ry'; failmap .* g2];
%failpoints(:,~any(failpoints))=[];
%plot3(passpoints(1,:), passpoints(2,:), passpoints(3,:), '.b')
hold on
%plot3(failpoints(1,:), failpoints(2,:), failpoints(3,:), '.r')
%For 2d plotting:
%s = surf(rx, ry, reshape(g2, size(rx)));
%s.EdgeColor = 'none';

[x, u, y, v, psi, r, delta, Fx] = decodeColocationVector(z);
plot(x,y, 'ob')
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), '--k')
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
%plot(sol1(:,1), sol1(:,3), 'b')
xlim([xlim1,xlim2])
ylim([ylim1,ylim2])
%view(0,90)

function [Lnormals, Rnormals] = calcNormals(leftBound, rightBound)
    n = size(leftBound, 2); % 246
    dLeft = leftBound(:,2:n) - leftBound(:,1:n-1);
    dRight = rightBound(:,2:n) - rightBound(:,1:n-1);
    Lnormals = [0, -1; 1, 0] * (dLeft ./ vecnorm(dLeft));
    Lnormals(:,n) = Lnormals(:,n - 1);
    Rnormals = [0, 1; -1, 0] * (dRight ./ vecnorm(dRight));
    Rnormals(:,n) = Rnormals(:,n - 1);
end

function g = inequal(pos, leftBound, rightBound, leftNormals, rightNormals)
    Idleft = knnsearch(leftBound', pos');
	Idright = knnsearch(rightBound', pos');

    rightCorr = rightBound(:, Idright');
    leftCorr = leftBound(:, Idleft');

    rightVectors = pos - rightCorr;
    leftVectors = pos - leftCorr;
    rightNormalsCorr = rightNormals(:, Idright');
    leftNormalsCorr = leftNormals(:, Idleft');

    g = -(dot(rightVectors, rightNormalsCorr, 1) .* dot(leftVectors, leftNormalsCorr, 1));

    %rightDirVectors = rightBound(:, 2:n) - rightBound(:, 1:n-1);
    %leftDirVectors = leftBound(:, 2:n) - leftBound(:, 1:n-1);
end