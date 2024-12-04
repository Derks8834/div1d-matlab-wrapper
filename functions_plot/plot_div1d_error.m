function plot_div1d_error(X, error,varargin)
% plot_div1d_error(x, error varargin)
% 
% INPUTS:
%  - X the DIV1D grid
%  - error (struct output from process_div1d_output.m)
%
% OUTPUT: produces figure

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

D.FontSize = 11;%2;
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

figure(P.fignum+3)
set(gcf, 'OuterPosition',[600, 100, 500, 500]);
subplot(2,2,1, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% plot the relative error
plot(X,error.density, 'LineWidth', P.LineWidth);
ylabel('rel err dens. eq.')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;
subplot(2,2,2, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% plot the relative pressure error
plot(X,error.momentum, 'LineWidth', P.LineWidth)
ylabel('rel err mom eq.')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;
title('stationarity and accuracy');

subplot(2,2,3, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% plot the relative qpar error
plot(X,error.energy, 'LineWidth', P.LineWidth)
ylabel('rel err ene eq.')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;
subplot(2,2,4, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% plot the relative neutral flux error
plot(X,error.neutral, 'LineWidth', P.LineWidth)
ylabel('rel error neutral eq.')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;

end