function fT = Jouanno(T, T_L, T_H,T_opt,steepness)
    if T < T_L || T > T_H
        fT = 0;
    else
        if T <= T_opt
            T_X = T_L;
        else
            T_X = T_H;
        end
        fT = exp(-steepness * ((T - T_opt) / (T_X - T))^2);
    end
end