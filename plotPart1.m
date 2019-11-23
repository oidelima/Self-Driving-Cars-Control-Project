%generatePart1

load TestTrack.mat

%sol1 = forwardIntegrateControlInput(ROB535ControlsProjectpart1input);

% Testing inequality constraint fitter
rng(0,'twister');
xlim1 = 200;
xlim2 = 1500;
rx = (xlim2-xlim1).*rand(10000,1) + xlim1;
ylim1 = -200;
ylim2 = 825;
ry = (ylim2-ylim1).*rand(10000,1) + ylim1;
[leftNormals, rightNormals] = calcNormals(TestTrack.bl, TestTrack.br);
% b = fitInequalConstraintToPoints(TestTrack.bl, 3)
% g1 = 200^2 - (rx - 900).^2 - (ry - 300).^2;
g2 = inequalityConstraintTrack([rx';ry'], TestTrack.bl, TestTrack.br);
passmap = g2 < 0;
failmap = g2 > 0;
passpoints = [passmap .* rx'; passmap .* ry'; passmap .* g2];
%passpoints(:,~any(passpoints))=[];
failpoints = [failmap .* rx'; failmap .* ry'; failmap .* g2];
%failpoints(:,~any(failpoints))=[];
plot3(passpoints(1,:), passpoints(2,:), passpoints(3,:), '.b')
hold on
plot3(failpoints(1,:), failpoints(2,:), failpoints(3,:), '.r')
%For 2d plotting:

plot3(TestTrack.cline(1,:), TestTrack.cline(2,:), zeros(1,246), '--k')
plot3(TestTrack.bl(1,:), TestTrack.bl(2,:), zeros(1,246), '.k')
plot3(TestTrack.br(1,:), TestTrack.br(2,:), zeros(1,246), '.k')
%plot(sol1(:,1), sol1(:,3), 'b')
xlim([200,1500])
ylim([-200, 825])
view(0,90)

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