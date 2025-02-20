function rate = eval_2D_Tne(coeff, T, ne)
    o = 0;
    for i = 0:8
        for j=0:8
            o = o + coeff(i+1, j+1)*(log(T).^i).*(log(ne*1e-14).^j);
        end
    end
    rate = 1e-6*exp(o);
end