dis = zeros(length(TestTrack.bl(1,:)),1);

for i = 1: length(length(TestTrack.bl(1,:)))
    dis(i) = sqrt((TestTrack.bl(1,i) - TestTrack.br(1,i))^2 + (TestTrack.bl(2,i) - TestTrack.br(2,i))^2);
end

function g = inequalityConstraintTrack(pos, leftBound, rightBound)
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
end

function [Lnormals, Rnormals] = calcNormals(pos,TestTrack.bl, TestTrack.br)
    n = size(leftBound, 2); % 246
    dLeft = leftBound(:,2:n) - leftBound(:,1:n-1);
    dRight = rightBound(:,2:n) - rightBound(:,1:n-1);

end