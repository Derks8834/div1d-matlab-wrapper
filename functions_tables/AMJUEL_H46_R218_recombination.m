function [recombination_AMJUEL] = AMJUEL_H46_R218_recombination(density,temperature)
% [recombination_AMJUEL] = recombination_AMJUEL_H46_R218(density,temperature)
% Inputs:
%   density in 1/m3  (may be vectors)
%   temperature in eV (may be vectors)
% Output:
%   recombination_AMJUEL in m3/s
%% effective recombination rate
% source AMJUEL (Note that densities in AMJUEL are in cm^-3)
% 4.6 reaction 2.1.8 
% G.l.Derks@differ.nl

recom_coef = [-2.858858570847E+01    2.068671746773E-02   -7.868331504755E-03  3.843362133859E-03   -7.411492158905E-04    9.273687892997E-05  -7.063529824805E-06    3.026539277057E-07   -5.373940838104E-09;
              -7.676413320499E-01    1.278006032590E-02   -1.870326896978E-02  3.828555048890E-03   -3.627770385335E-04    4.401007253801E-07   1.932701779173E-06   -1.176872895577E-07    2.215851843121E-09;
               2.823851790251E-03   -1.907812518731E-03    1.121251125171E-02 -3.711328186517E-03    6.617485083301E-04   -6.860774445002E-05   4.508046989099E-06   -1.723423509284E-07    2.805361431741E-09;
              -1.062884273731E-02   -1.010719783828E-02    4.208412930611E-03 -1.005744410540E-03    1.013652422369E-04   -2.044691594727E-06  -4.431181498017E-07    3.457903389784E-08   -7.374639775683E-10;
               1.582701550903E-03    2.794099401979E-03   -2.024796037098E-03  6.250304936976E-04   -9.224891301052E-05    7.546853961575E-06  -3.682709551169E-07    1.035928615391E-08   -1.325312585168E-10;
              -1.938012790522E-04    2.148453735781E-04    3.393285358049E-05 -3.746423753955E-05    7.509176112468E-06   -8.688365258514E-07   7.144767938783E-08   -3.367897014044E-09    6.250111099227E-11;
               6.041794354114E-06   -1.421502819671E-04    6.143879076080E-05 -1.232549226121E-05    1.394562183496E-06   -6.434833988001E-08  -2.746804724917E-09    3.564291012995E-10   -8.551708197610E-12;
               1.742316850715E-06    1.595051038326E-05   -7.858419208668E-06  1.774935420144E-06   -2.187584251561E-07    1.327090702659E-08  -1.386720240985E-10   -1.946206688519E-11    5.745422385081E-13;
              -1.384927774988E-07   -5.664673433879E-07    2.886857762387E-07 -6.591743182569E-08    8.008790343319E-09   -4.805837071646E-10   6.459706573699E-12    5.510729582791E-13   -1.680871303639E-14 ];
 
%
% selected densities 
% note normalization by 10^14 of densities in AMJUEL fits
%density = [ 1.0e+18 1.0e+20 1.0e+21 ]/1.0e14;
density = density/1.0e14; 

recombination_AMJUEL = 0.0 * density' * temperature;
for i = 1: length(temperature) % temperature loop
    for j = 1: length(density) % density loop
        for m = 1: 9
            for n = 1: 9
                recombination_AMJUEL(j,i) = recombination_AMJUEL(j,i) + recom_coef(n,m) * log(density(j))^(m-1) * log(temperature(i))^(n-1);
            end
        end
    end
end
recombination_AMJUEL = exp(recombination_AMJUEL) * 1.0d-6;

end