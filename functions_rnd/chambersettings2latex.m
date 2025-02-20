function chambersettings2latex(input)
% function to plot latex table of DIV1D reservoir inputs with core-SOL
% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% Jan 2025
%
% 

% shots and puff should be present to give the list 
% (otherwise the settings make no sense)
shots = input.shots; 
puff = input.puff; 
F = fields(input);
fprintf('|c|')
for ii = 1:length(shots)
fprintf('r|')
end
fprintf('c| \n')

for ii = 1:length(F)

    tmpdat = input.(F{ii});
    if ~strcmp(F{ii},'shots')
    [n,m] = size(tmpdat);
    % depending on dimensions print matrix or just row
    if n > 1
        for jj = 1:m
        fprintf('\\verb| %s %s| ',F{ii},num2str(jj))
        fprintf(' & %1.2e',tmpdat(:,jj))
        fprintf('& [-]\\\\ \\hline \n')
        end
    else
    fprintf('\\verb| %s | ',F{ii})
    fprintf(' & %1.2e ',tmpdat)
    fprintf('& [-]\\\\ \\hline \n')
    end
    else
    fprintf('\\verb| %s | ',F{ii})
    fprintf(' & %i ',tmpdat)
    fprintf('& [-]\\\\ \\hline \n')
    end
    
end



fprintf('\n \n')

end