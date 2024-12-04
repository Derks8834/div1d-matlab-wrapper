function [elastic_AMJUEL] = AMJUEL_P38_R21_tot_elast_rate(temperature)
%AMJUEL_P38_R21_TOT_ELAST_RATE Summary of this function goes here
% source AMJUEL page 38 2.1 (total elastic rate coefficient)
%
% INPUT: temperature in eV.
% Output: elastic_AMJEL m3/(s)
% G.l.Derks@differ.nl


% coefficients of fit:
b_coef = [ -1.833882000000E+01 2.368705000000E-01 -1.469575000000E-02 -1.139850000000E-02 6.379644000000E-04 3.162724000000E-04 -6.681994000000E-05 3.812123000000E-06 8.652321000000E-09 ];
for i = 1: length(temperature)
    for j = 1: 9
        elastic_AMJUEL(i) = elastic_AMJUEL(i) + b_coef(j)*log(temperature(i))^(j-1);
    end
end
elastic_AMJUEL = exp(elastic_AMJUEL)*1.0d-6;

end

