function b = fitInequalConstraintToPoints(boundary1, centerline, boundary2, order)
    syms X Y
    depth = -1;
    x = [boundary1(1,:)'; centerline(1,:)'; boundary2(1,:)'];
    y = [boundary1(2,:)'; centerline(2,:)'; boundary2(2,:)'];
    m1 = size(boundary1, 2);
    m2 = size(centerline, 2);
    m3 = size(boundary2, 2);
    m = m1 + m2 + m3;
    A1 = ones(m, 1);
    A2 = zeros(m, order);
    A3 = zeros(m, order);
    for i=1:order
        A2(:,i) = x.^i;
        A3(:,i) = y.^i;
        S(1+i) = X^i;
        S(1+i+order) = Y^i;
    end
    A = [A1,A2,A3];
    S(1) = 1;
    levels = [zeros(m1,1);depth*ones(m2,1);zeros(m3,1)];
    x
    y
    levels
    A
    b = A\levels;
    A*b
end