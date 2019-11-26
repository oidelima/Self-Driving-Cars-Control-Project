function z = encodeColocationVector(x, u, y, v, psi, r, delta, Fx)
    n_minus1 = size(delta, 1);
    l1 = (n_minus1+1)*6;
    l2 = n_minus1*2;
    z1 = reshape([x, u, y, v, psi, r]',[l1,1]);
    z2 = reshape([delta, Fx]',[l2,1]);
    % Final size of z is [l1 + l2, 1] = [6n + 6 + 2n, 1] = [8n + 6, 1]
    z = [z1; z2];
end