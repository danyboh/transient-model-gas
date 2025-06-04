function [DRo, Rom, z, A1, A2, A3, A4] = fdens(p, T, Rom)
    global Ropk Tpk c R
    r = 10;
    S = [7 6 6 5 5 4 3 3 3 2];
    Rop = Rom / Ropk;
    Tp = T / Tpk;

    z = 1; A1 = 0; A2 = 0; A3 = 0; A4 = 0;

    for k = 1:r
        l_range = 0:S(k);
        cl = c(k, l_range + 1);
        Tp_pow = Tp .^ l_range;
        Rop_pow = Rop ^ k;

        z = z + sum(cl .* Rop_pow ./ Tp_pow);
        A1 = A1 + sum((k + 1) * cl .* Rop_pow ./ Tp_pow);
        A2 = A2 - sum((l_range - 1) .* cl .* Rop_pow ./ Tp_pow);
        A3 = A3 - sum(l_range .* (l_range - 1) / k .* cl .* Rop_pow ./ Tp_pow);
        A4 = A4 - sum(l_range .* cl .* Rop_pow ./ (Tp_pow .* Tp * Tpk));
    end

    DRo = (1e3 * p - R * T * z * Rom) / (R * T * (1 + A1));
    Rom = Rom + DRo;
end
