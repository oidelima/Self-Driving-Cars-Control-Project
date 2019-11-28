function [z, startpts] = generateMultiplePart1(nsteps, niter, subdivision_num, gradient_precision)
    global n dt
    z0 = repmat([287, 5, -176, 0, 2, 0], n, 1);
    u0 = repmat([0 0], n-1, 1);

    x_total = [];
    u_total = [];
    y_total = [];
    v_total = [];
    psi_total = [];
    r_total = [];
    delta_total = [];
    Fx_total = [];

    for iteration = 1:niter
        startpts(iteration, :) = [z0(1,1), z0(1,2), z0(1,3), z0(1,4), z0(1,5), z0(1,6)];
        x0 = encodeColocationVector(z0(:,1), z0(:,2), z0(:,3), z0(:,4), z0(:,5), z0(:,6), u0(:,1), u0(:,2));
        iteration
        tic
        z_out = generatePart1(nsteps, x0, subdivision_num, gradient_precision);
        toc
    
        [x, u, y, v, psi, r, delta, Fx] = decodeColocationVector(z_out);
        z0 = repmat([x(end,1), u(end,1), y(end,1), v(end,1), psi(end,1), r(end,1)], n, 1);
        u0 = repmat([delta(end,1), Fx(end,1)], n-1, 1);

        if iteration==niter
            x_total = [x_total;x(1:end)];
            u_total = [u_total;u(1:end)];
            y_total = [y_total;y(1:end)];
            v_total = [v_total;v(1:end)];
            psi_total = [psi_total;psi(1:end)];
            r_total = [r_total;r(1:end)];
            delta_total = [delta_total;delta];
            Fx_total = [Fx_total;Fx];
        else
            x_total = [x_total;x(1:end-1)];
            u_total = [u_total;u(1:end-1)];
            y_total = [y_total;y(1:end-1)];
            v_total = [v_total;v(1:end-1)];
            psi_total = [psi_total;psi(1:end-1)];
            r_total = [r_total;r(1:end-1)];
            delta_total = [delta_total;delta];
            Fx_total = [Fx_total;Fx];
        end
    end
    z = encodeColocationVector(x_total, u_total, y_total, v_total, psi_total, r_total, delta_total, Fx_total);
end