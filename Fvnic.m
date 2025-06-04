function [Kvnic, z, A1, A2, A3, A4, Mm, Rom_r] = Fvnic(p, T, x)
    persistent cached_static;
    global O Rok M lambda ksi Tk 
    if isempty(cached_static)
        cached_static = struct();
        dat_vnic;
        cached_static.R = 8.31451;
        cached_static.Tpk_ready = false;
    end
    R = cached_static.R;

    if length(x) == 20
        
        x_tmp = zeros(1,8);
        x_tmp(1:5) = x(1:5);
        x_tmp(6) = sum(x(6:10)) + sum(x(14:18));
        x_tmp(7) = x(11);
        x_tmp(8) = x(15);
        x = x_tmp;
    elseif length(x) == 18
        x_tmp = x(1:8);
        x_tmp(4) = x_tmp(4) + sum(x(9:13));
        x_tmp(6) = x_tmp(6) + sum(x(14:18));
        x = x_tmp;
    end

    if ~cached_static.Tpk_ready
        global O Rok M lambda ksi Oij Vk Tkij
        n = 8;
        Oij = zeros(n); Vk = zeros(n); Tkij = zeros(n);
        for i = 1:n
            for j = 1:n
                Oij(i,j) = (O(i)*M(i)/Rok(i) + O(j)*M(j)/Rok(j)) / (M(i)/Rok(i) + M(j)/Rok(j));
                Vk(i,j) = (1 - lambda(i,j)) * (((M(i)/Rok(i))^(1/3) + (M(j)/Rok(j))^(1/3)) / 2)^3;
                Tkij(i,j) = (1 - ksi(i,j)) * sqrt(Tk(i)*Tk(j));
            end
        end
        Xi = repmat(x,8,1); Xj = Xi';
        cached_static.Ropk = 1 / sum(sum(Xi .* Xj .* Vk));
        cached_static.Om = cached_static.Ropk * sum(sum(Xi .* Xj .* Vk .* Oij));
        cached_static.Tkm = sum(sum(Xi .* Xj .* Vk .* Tkij.^2));
        cached_static.Tpk = sqrt(cached_static.Tkm * cached_static.Ropk);
        cached_static.Tpk_ready = true;
    end

    Ropk = cached_static.Ropk; Tpk = cached_static.Tpk; Om = cached_static.Om;
    Ppk = 1e-3 * R * Ropk * Tpk * (0.28707 - 0.05559 * Om);

    Pp = p / Ppk;
    Rom = 9e3 * p / (R * T * (1.1 * Pp + 0.7));

    [DRo, Rom, z, A1, A2, A3, A4] = fdens(p, T, Rom);
    iter = 1; max_iter = 20; tol = 1e-4;
    while abs(DRo / Rom) >= tol && iter < max_iter
        [DRo, Rom, z, A1, A2, A3, A4] = fdens(p, T, Rom);
        iter = iter + 1;
    end
    if iter >= max_iter
        warning('Fvnic not converged: p = %.2f, T = %.2f', p, T);
    end

    Rom_r = Rom;
    Pc = 0.101325; Tc = 293.15;
    Pp = Pc / Ppk;
    Rom = 9e3 * Pc / (R * Tc * (1.1 * Pp + 0.7));
    [DRo, Rom, zc] = fdens(Pc, Tc, Rom);
    iter = 1;
    while abs(DRo / Rom) >= tol && iter < max_iter
        [DRo, Rom, zc] = fdens(Pc, Tc, Rom);
        iter = iter + 1;
    end

    Kvnic = z / zc;
    Mm = sum(x .* M);
end
