function [x, u, y, v, psi, r, delta, Fx] = decodeColocationVector(z)
    % 8n + 6 = size(z,1)
    n_minus1 = (size(z, 1) - 6) / 8;
    l1 = (n_minus1+1)*6;
    l2 = n_minus1*2;
    z = reshape(z, [l1+l2,1]);
    x = z(1:6:l1);
    u = z(2:6:l1);
    y = z(3:6:l1);
    v = z(4:6:l1);
    psi = z(5:6:l1);
    r = z(6:6:l1);
    delta = z(l1+1:2:l1+l2);
    Fx = z(l1+2:2:l1+l2);
end