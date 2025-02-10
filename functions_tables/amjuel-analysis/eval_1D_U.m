function U = eval_1D_U(coeffs, r)
%  Inputs: AMJUEL coefficients
%             r: distance in Bohr radii
%  Output: U: potential energy in eV

    epsilon = coeffs(1);
    g1 = coeffs(2);
    g2 = coeffs(3); 
    rm = coeffs(4); 
    
    rho = r./rm;
    g = zeros(1,length(rho));
    g(rho<1) = g1;
    g(rho>=1) = g1*g2;

    U = epsilon * (exp(2*g.*(1-rho))-2*exp(g.*(1-rho)));
end
