function input2latex(input)
% simple function to plot a latex table for DIV1D inputs
% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

shots = input.shots;
puff = input.puff;
L = input.L;
q_par_X = input.q_par_X;
n_x = input.n_x;
cratio = input.cratio;
sintheta = input.sintheta;

fprintf('|c|r|r|r|r|r|r|r|r|r|c| \n')
fprintf('SOLPS')
fprintf(' & %i ',shots)
fprintf('& [-]\\\\ \\hline \n')

fprintf('puff ')
fprintf(' & %i ',puff/10^20)
fprintf('& [**20 $D$/s] \\\\ \n')

fprintf('\\verb|L|')
fprintf('& %2.2f ',L)
fprintf('& [m] \\\\ \n')

fprintf( '\\verb|q_parX|')
fprintf('& %2.2f ',q_par_X/10^6)
fprintf('& [MW/m2] \\\\ \n')

fprintf( '\\verb|n_X|')
fprintf('& %2.2f ',n_x/10^19)
fprintf('& [**19/m3] \\\\ \n')

fprintf( '\\verb|car_con|')
fprintf('& %1.3f ',cratio)
fprintf('& [-] \\\\ \n')

fprintf('\\verb|sintheta|')
fprintf('& %1.3f ',sintheta)
fprintf('& [-] \\\\ \\hline \n')
try
% radlos = input.radlos;
% fprintf('\\verb|radloss|')
% fprintf('& %1.3f ',radlos)
% fprintf('& [-] \\\\  \n')    

% rad_exp_los = input.rad_exp_los;
% fprintf('\\verb|rad_exp_los|')
% fprintf('& %1.3f ',rad_exp_los)
% fprintf('& [-] \\\\  \n')  

Ru_Rt_mwid = input.Ru_Rt_mwid;
fprintf('\\verb|Ru_Rt_mwid|')
fprintf('& %1.3f ',Ru_Rt_mwid)
fprintf('& [-] \\\\ \\hline \n')  

flux_exp = input.flx_exp;
fprintf('\\verb|ef|')
fprintf('& %1.3f ',flux_exp)
fprintf('& [-] \\\\ \\hline \n')  

n_b = input.n_b;
fprintf('\\verb|n\_b|')
fprintf('& %1.3f ',n_b)
fprintf('& [-] \\\\ \\hline \n')  


catch
end

fprintf('\n \n')

end