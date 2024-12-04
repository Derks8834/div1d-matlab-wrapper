function plot_div1d_error_v600(X, error,varargin)
% plot_div1d_error(x, error varargin)
% 
% INPUTS:
%  - X the DIV1D grid
%  - error (struct output from process_div1d_output.m)
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
errorfields = fields(error);
for ie = 1:length(errorfields)
    figure(ie+ P.fignum);
    efields = fields(error.(errorfields{ie}));
    N = length(efields);
    for ip = 1:N
        subplot(N,1,ip)
        plot(X,error.(errorfields{ie}).(efields{ip}));
        ylabel(strrep(efields{ip},'_',' '));
        if P.hold ==1; hold on; end
        if P.flip ==1; set(gca,'XDir','reverse'); end     
    end
    subplot(N,1,1)
    title(strcat('errors: ',errorfields{ie} ));
end
