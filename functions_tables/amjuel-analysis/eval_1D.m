function rate = eval_1D(coeff, T)
    o = 0;
    for i = 0:8
        o = o + coeff(i+1).*(log(T).^i);
    end
    rate = 1e-6*exp(o);
end
