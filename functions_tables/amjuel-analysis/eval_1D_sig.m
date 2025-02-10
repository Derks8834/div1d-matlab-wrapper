function rate = eval_1D_sig(coeff, E)
    o = 0;
    for i = 0:8
        o = o + coeff(i+1).*(log(E).^i);
    end
    rate = exp(o);
end