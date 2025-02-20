function [] = plot_div1d_bal_v600(X,bal,varargin)
% function -> plot_div1d_bal_v600(X,bal,varargin) 
% plots the terms in the balance equations of div1d 
% terms are calculated by process_div1d_output
% Input: X is the xgrid of div1d
% Input: bal struct containing balance for density, energy, etc.
%
% OUTPUT: produces figure

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2024
D.FontSize = 11;
D.fignum = 16;
D.LineWidth = 1;
D.hold = 0;
D.avg = 0;
D.flip = 1;

P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end
if P.flip ==1
    X = max(X) - X;
end  

figure(P.fignum);
balfields = fields(bal);
N = length(balfields);
for ie = 1:length(balfields)
    balent = fields(bal.(balfields{ie}));
    ylimits = [0  0.01];
    subplot(1,N,ie)
    for ip = 1:length(balent)  
         ytmp = bal.(balfields{ie}).(balent{ip});
         nx = length(ytmp);
        plot(X,ytmp); hold on;
           if P.hold ==1; hold on; end
        if P.flip ==1; set(gca,'XDir','reverse'); end     
        ylimits(1) = min(ylimits(1),min(ytmp(3:nx-2)));
        ylimits(2) = max(ylimits(2),max(ytmp(3:nx-2)));
    end

    ylim(ylimits);
    ylabel(balfields{ie})
    legend(balent{:});
    hold off;
    if P.flip ==1; xlabel('distance to target [m]'); end

end
if P.flip ==1; xlabel('distance to target [m]'); 
else
xlabel('distance to upstream [m')
end    

end

