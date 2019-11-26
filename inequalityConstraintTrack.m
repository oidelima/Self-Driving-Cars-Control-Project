function g = inequalityConstraintTrack(pos, leftBound, rightBound)
    upperbound = 10000;

    [leftNormals, rightNormals] = calcNormals(leftBound, rightBound);
    Idleft = knnsearch(leftBound', pos');
	Idright = knnsearch(rightBound', pos');

    rightCorr = rightBound(:, Idright');
    leftCorr = leftBound(:, Idleft');

    rightVectors = pos - rightCorr;
    leftVectors = pos - leftCorr;
    rightNormalsCorr = rightNormals(:, Idright');
    leftNormalsCorr = leftNormals(:, Idleft');
    

    g = -(dot(rightVectors, rightNormalsCorr, 1) .* dot(leftVectors, leftNormalsCorr, 1));
    %g = g_temp .* (g_temp < upperbound) + upperbound * (g_temp >= upperbound);
end

function [Lnormals, Rnormals] = calcNormals(leftBound, rightBound)
    n = size(leftBound, 2); % 246
    dLeft = leftBound(:,2:n) - leftBound(:,1:n-1);
    dRight = rightBound(:,2:n) - rightBound(:,1:n-1);
    Lnormals = [0, -1; 1, 0] * (dLeft ./ vecnorm(dLeft));
    Lnormals(:,n) = Lnormals(:,n - 1);
    Rnormals = [0, 1; -1, 0] * (dRight ./ vecnorm(dRight));
    Rnormals(:,n) = Rnormals(:,n - 1);
end