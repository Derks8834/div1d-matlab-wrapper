function input2latex_stag(input)
% function to plot latex table of DIV1D inputs with core-SOL
% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

shots = input.shots;
%puff = input.puff;
L = input.L;
L_core_SOL = input.L_core_SOL;
Q_core = input.Q_core;
Gamma_core = input.Gamma_core;
sintheta_mid = input.sintheta_mid;
sintheta_tar = input.sintheta_tar;
lambda_mid = input.width_mid;
lambda_tar = input.width_tar;
maj_rad= input.major_radius;
expansion_total = input.expansion_total3;
expansion_flux = input.expansion_flux;
expansion_trans=input.expansion_trans;
neutral_density_div_E = input.neutral_density_div_E;
neutral_density_cor_E = input.neutral_density_cor_E;
try neutral_density_div_I = input.neutral_density_div_I;
neutral_density_mdiv_I = input.neutral_density_mdiv_I;
neutral_density_mdiv_E = input.neutral_density_mdiv_E;
catch; end;
neutral_density_cor = input.neutral_density_cor;
neutral_density_div = input.neutral_density_div;
carbon_concentration = input.carbon_concentration;

try 
core2solflux = input.flux_i_core2sol;
sol2divflux = input.flux_i_sol2div;
mfp_ion_h4215 = input.mfp_ion_h4215; % mean free path ionization
gamma_omp = input.gamma_omp; % shielding factor
lambda_te = input.lambda_te;
catch

end


fprintf('|c|r|r|r|r|r|r|r|r|r|c| \n')
fprintf('SOLPS')
fprintf(' & %i ',shots)
fprintf('& [-]\\\\ \\hline \n')

fprintf('\\verb|L_{core}|')
fprintf('& %2.2f ',L_core_SOL)
fprintf('& [m] \\\\ \n')

fprintf('\\verb|L_{div}|')
fprintf('& %2.2f ',L-L_core_SOL)
fprintf('& [m] \\\\ \n')

fprintf('\\verb|major_radius|')
fprintf('& %1.3f ',maj_rad)
fprintf('& [m] \\\\ \\hline \n')

fprintf('\\verb|lambda_core|')
fprintf('& %1.3f ',lambda_mid)
fprintf('& [m] \\\\ \\hline \n')

fprintf('\\verb|lambda_div|')
fprintf('& %1.3f ',lambda_tar)
fprintf('& [m] \\\\ \\hline \n')

fprintf('\\verb|sintheta_core|')
fprintf('& %1.3f ',sintheta_mid)
fprintf('& [-] \\\\ \\hline \n')

fprintf('\\verb|sintheta_div|')
fprintf('& %1.3f ',sintheta_tar)
fprintf('& [-] \\\\ \\hline \n')

fprintf( '\\verb|Q_{core}|')
fprintf('& %2.2f ',Q_core/10^6)
fprintf('& [MW] \\\\ \n')

fprintf( '\\verb|Gamma_{core}|')
fprintf('& %2.2f ',Gamma_core/10^21)
fprintf('&[#10^{21}s^{-1}] \\\\ \n')


fprintf( '\\verb|car_con|')
fprintf('& %1.3f ',carbon_concentration)
fprintf('& [-] \\\\ \n')

fprintf( '\\verb|nb_core|')
fprintf('& %1.2e ',neutral_density_cor_E)
fprintf('& [m^{-3}] \\\\ \n')
fprintf( '\\verb|nb_div|')
fprintf('& %1.2e ',neutral_density_div_E)
fprintf('& [m^{-3}] \\\\ \n')

try
fprintf( '\\verb|nb_mdiv_I|')
fprintf('& %1.2e ',neutral_density_mdiv_I)
fprintf('& [m^{-3}] \\\\ \n')
fprintf( '\\verb|nb_div_I|')
fprintf('& %1.2e ',neutral_density_div_I)
fprintf('& [m^{-3}] \\\\ \n')

fprintf( '\\verb|nb_mcore|')
fprintf('& %1.2e ',neutral_density_mdiv_E)
fprintf('& [m^{-3}] \\\\ \n')
fprintf( '\\verb|nb_mdiv_E|')
fprintf('& %1.2e ',neutral_density_mcor_E)
fprintf('& [m^{-3}] \\\\ \n')
catch end

fprintf( '\\verb|nn_core|')
fprintf('& %1.2e ',neutral_density_cor)
fprintf('& [m^{-3}] \\\\ \n')
fprintf( '\\verb|nn_div|')
fprintf('& %1.2e ',neutral_density_div)
fprintf('& [m^{-3}] \\\\ \n')

fprintf( '\\verb|fl_exp|')
fprintf('& %1.2f ',expansion_flux)
fprintf('& [-] \\\\ \n')
fprintf( '\\verb|tr_exp|')
fprintf('& %1.2f ',expansion_trans)
fprintf('& [-] \\\\ \n')

fprintf( '\\verb|Gam_core2sol|')
fprintf('& %1.2e ',core2solflux)
fprintf('& [#/s] \\\\ \n')

fprintf( '\\verb|Gam_sol2div|')
fprintf('& %1.2e ',sol2divflux)
fprintf('& [#/s] \\\\ \n')

try
fprintf( '\\verb|lambda_te_omp|')
fprintf('& %1.2f ',lambda_te)
fprintf('& [m] \\\\ \n')
catch; end

try
fprintf( '\\verb|mfp_ion/lambda_te_omp|')
fprintf('& %1.2f ',gamma_omp)
fprintf('& [-] \\\\ \n')
catch; end

fprintf('\n \n')

end