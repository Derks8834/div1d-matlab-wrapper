function [charge_exchange_AMJUEL] = AMJUEL_H219_R318_charge_exchange(temperature,mass)
%% effective charge exchange rate from AMJUEL as function of temperature
% temperature in eV input.
% source AMJUEL (Note that densities in AMJUEL are in cm^-3)
% 2.19 reaction 3.1.8 

mass_h = 1.67262192369*10^(-27); % kg
temperature = max(temperature*mass_h/mass,0.1);
cx_coef = [-1.850280000000E+01  3.708409000000E-01  7.949876000000E-03 ...
           -6.143769000000E-04 -4.698969000000E-04 -4.096807000000E-04 ...
            1.440382000000E-04 -1.514243000000E-05  5.122435000000E-07];

cx_AMJUEL = 0.0 * temperature;
for i = 1: length(temperature) % temperature loop
    for n = 1: 9
        cx_AMJUEL(i) = cx_AMJUEL(i) + cx_coef(n) * log(temperature(i))^(n-1);
    end
end
charge_exchange_AMJUEL = exp(cx_AMJUEL) * 1.0d-6;
end