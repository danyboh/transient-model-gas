% === CpZ_cached.m ===
function [Cp_out, Z_out] = CpZ_cached(p, T, x_mol)
    persistent cache;
    if isempty(cache)
        cache = containers.Map;
    end

    key = sprintf('%.0f_%.0f', round(p), round(T));
    if isKey(cache, key)
        val = cache(key);
        Cp_out = val(1);
        Z_out = val(2);
    else
        [Cp_out, ~, ~] = Cp_Vnic(p, T, x_mol);
        [~, Z_out, ~, ~, ~, ~, ~] = Fvnic(p, T, x_mol);
        cache(key) = [Cp_out, Z_out];
    end
end
