function rate = eval_2D_TE(coeff, T, E)
    o = 0;
    for i = 0:8
        for j=0:8
            o = o + coeff(i+1, j+1).*(log(T).^i).*(log(E).^j);
        end
    end
    rate = 1e-6*exp(o);
end