function plotdiv1d_v600(o,i,varargin)
% plot_div1d_output(div1doutput,div1dinput,proc,varargin)
% div1doutput and input and proc obtained by proces_div1d_output.m
% 
% varargin
% D.solps = 0;
% D.S2PM = 0;
% D.D2PM = 0;
% D.cratio = 0; % in percentage
% D.dxmin = 0.05; 
% D.Ntime = length(div1doutput.time);
% D.flip = 1;
% D.FontSize = 11;%2;
% D.fignum = 16;
% D.LineWidth = 1;
% D.hold = 0;
% D.avg = 1;
% D.save = 0;
% D.generalposition = [100 100 200 250];
% D.title = '' ;
% D.Dtitle = '' ;

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

D.solps = 0;
D.S2PM = 0;
D.D2PM = 0;
D.cratio = 0; % in percentage
D.dxmin = 0.05; 
D.Ntime = length(o.time); %,i.numerics.ntime);
D.flip = 1;
D.FontSize = 11;%2;
D.fignum = 10;
D.LineWidth = 1;
D.hold = 1;
D.avg = 1;
D.save = 0;
D.generalposition = [100 100 200 250];
D.title = '' ;
D.Dtitle = '' ;

P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end

e_charge = 1.6021766341e-19;
%B_field = input.flx_exp_prf;
Nt = P.Ntime;
Nx = i.numerics.nx;
dxmin = i.numerics.dxmin;
x = i.grid.x;
xcb = i.grid.xcb;
% momentum % post processing here? 
%P_kin = 2.0*o.density(Nt,:).*o.temperature(Nt,:)*e_charge;
%P_tot =  P_kin + o.Gamma_mom(Nt,1:end-1); 

%% General
cmap = lines;
xstr = 'distance from upstream [m]';
if P.flip ==1
    x = max(x) - x;
    xcb = max(xcb) - xcb;
    xstr = 'distance to target [m]';
end

%% Profiles
figure(P.fignum)
% temperature
subplot(2,3,1, 'FontSize', P.FontSize)
plot(x,o.temperature(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('temperature [eV]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% density
subplot(2,3,2, 'FontSize', P.FontSize)
semilogy(x,o.density(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('density [/m3]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% velocity
subplot(2,3,3, 'FontSize', P.FontSize)
plot(x,o.velocity(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('velocity [m/s]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end
% atoms
subplot(2,3,4, 'FontSize', P.FontSize)
semilogy(x,o.neutral_density(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('atoms [/m3]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% atom velocity
subplot(2,3,5, 'FontSize', P.FontSize)
plot(x,o.neutral_velocity(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('atom velocity [m/s]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% molecule
subplot(2,3,6, 'FontSize', P.FontSize)
try
semilogy(x,o.molecule(Nt,:), 'LineWidth', P.LineWidth);  grid on;
catch
semilogy(x,o.molecule_density(Nt,:), 'LineWidth', P.LineWidth);  grid on;
end
if P.hold == 1;  hold on; end
ylabel('molecule [m/s]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

set(gcf, 'OuterPosition',[10, 10, 700, 500]);


%% FLUXES
figure(P.fignum +1)
%set(gcf, 'OuterPosition',[100, 600, 500, 500]);
% q parallel
subplot(2,3,1, 'FontSize', P.FontSize)
plot(xcb,o.q_parallel(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('q [W/m2]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% density flux
subplot(2,3,2, 'FontSize', P.FontSize)
plot(xcb,o.Gamma_n(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on;  end
ylabel('Gamman n []');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% velocity momentum flux
subplot(2,3,3, 'FontSize', P.FontSize)
plot(xcb,o.Gamma_mom(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('Gamma mom []');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end
% atom flux
subplot(2,3,4, 'FontSize', P.FontSize)
try
plot(xcb,o.Gamma_neutral(Nt,:), 'LineWidth', P.LineWidth);  grid on;
catch
    plot(xcb,o.neutral_flux(Nt,:), 'LineWidth', P.LineWidth);  grid on;
end

if P.hold == 1;  hold on; end
ylabel('Gamma atom []');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% atom momentum flux
subplot(2,3,5, 'FontSize', P.FontSize)
try
plot(xcb,o.Gamma_mom_neutral(Nt,:), 'LineWidth', P.LineWidth);  grid on;
catch
plot(xcb,o.Gamma_mom_n(Nt,:), 'LineWidth', P.LineWidth);  grid on;
end
if P.hold == 1;  hold on; end
ylabel('gamma atom velocity []');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% molecule flux
subplot(2,3,6, 'FontSize', P.FontSize)
try
plot(xcb,o.Gamma_molecule(Nt,:), 'LineWidth', P.LineWidth);  grid on;
catch
plot(xcb,o.molecule_flux(Nt,:), 'LineWidth', P.LineWidth);  grid on;
end
if P.hold == 1;  hold on; end
ylabel('Gamma molecule [/m2s]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

set(gcf, 'OuterPosition',[10, 10, 700, 500]);


%% SOURCES
figure(P.fignum +2)
%set(gcf, 'OuterPosition',[600, 600, 500, 500]);
% temperature
subplot(2,3,1, 'FontSize', P.FontSize)
plot(x,o.Source_Q(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('Q []');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% particles
subplot(2,3,2, 'FontSize', P.FontSize)
plot(x,o.Source_n(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('Source n e []');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% velocity
subplot(2,3,3, 'FontSize', P.FontSize)
plot(x,o.Source_v(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('Source v []');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end
% atoms
subplot(2,3,4, 'FontSize', P.FontSize)
plot(x,o.Source_neutral(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('Source atoms []');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% atom velocity
subplot(2,3,5, 'FontSize', P.FontSize)
plot(x,o.Source_vn(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('Source v atom []');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

% molecule
subplot(2,3,6, 'FontSize', P.FontSize)
plot(x,o.Source_molecule(Nt,:), 'LineWidth', P.LineWidth);  grid on;
if P.hold == 1;  hold on; end
ylabel('Source molecule [m/s]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

set(gcf, 'OuterPosition',[10, 10, 700, 500]);


%% plot cross-sol fluxes
figure(P.fignum +3) 
subplot(221, 'FontSize', P.FontSize)
try
plot(x,o.core2sol_flux(Nt,:));
catch
disp('core2sol was not written in old branch')
end
if P.hold == 1;  hold on; end
ylabel('cor2sol flux [m2/s]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end
subplot(222, 'FontSize', P.FontSize)
plot(x,o.core2sol_mol(Nt,:));
if P.hold == 1;  hold on; end
ylabel('core2sol mol [m2/s]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end
subplot(223, 'FontSize', P.FontSize)
plot(x,o.extern2sol_flux(Nt,:));
if P.hold == 1;  hold on; end
ylabel('ext2sol flux [m2/s]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end
subplot(224, 'FontSize', P.FontSize)
plot(x,o.extern2sol_mol(Nt,:));
if P.hold == 1;  hold on; end
ylabel('ext2sol mol [m2/s]');
xlabel(xstr); 
if P.flip ==1; set(gca,'XDir','reverse');end

set(gcf, 'OuterPosition',[10, 10, 500, 500]);

%% plot div1d reservoirs
figure(P.fignum+4);
subplot(331)% DFR 1
plot(o.time,o.extern_neutral_density(:,2));
if P.hold == 1;  hold on; end
title('inner dfr')
subplot(332) % CFR
plot(o.time,o.extern_neutral_density(:,3));
if P.hold == 1;  hold on; end
title('ATOMS cfr')
subplot(333) % DFR2
plot(o.time,o.extern_neutral_density(:,4));
if P.hold == 1;  hold on; end
title('outer dfr')
subplot(334) % inner target sol
plot(o.time,mean(o.neutral_density(:, 1: max(i.grid.i_xpoint(1),1) ),2) )
if P.hold == 1;  hold on; end
title('inner sol')
subplot(335) % core sol
plot(o.time, mean(o.neutral_density(:, max(i.grid.i_xpoint(1),1) : max(i.grid.i_xpoint(2),1) ),2) )
if P.hold == 1;  hold on; end
title('core sol')
subplot(336) % outer target sol
plot(o.time, mean(o.neutral_density(:, max(i.grid.i_xpoint(2),1) : end ),2) )
if P.hold == 1;  hold on; end
title('outer sol')
subplot(337) % in PFR 
plot(o.time,o.extern_neutral_density(:,1));
if P.hold == 1;  hold on; end
title('inner pfr')
xlabel('time (s)')
subplot(338) % core
plot(o.time,o.core_neutral_density)
if P.hold == 1;  hold on; end
title('core atom')
xlabel('time (s)')
subplot(339) % out PFR
plot(o.time,o.extern_neutral_density(:,5));
if P.hold == 1;  hold on; end
title('outer pfr')
xlabel('time (s)')

figure(P.fignum+5);
subplot(331)% DFR 1
plot(o.time,o.extern_molecule_density(:,2));
if P.hold == 1;  hold on; end
title('inner dfr')
subplot(332) % CFR
plot(o.time,o.extern_molecule_density(:,3));
if P.hold == 1;  hold on; end
title('MOLECULES cfr')
subplot(333) % DFR2
plot(o.time,o.extern_molecule_density(:,4));
if P.hold == 1;  hold on; end
title('outer dfr')
subplot(334) % inner target sol
plot(o.time,mean(o.molecule(:, 1: max(i.grid.i_xpoint(1),1) ),2) )
if P.hold == 1;  hold on; end
title('inner sol')
subplot(335) % core sol
plot(o.time, mean(o.molecule(:, max(i.grid.i_xpoint(1),1) : max(i.grid.i_xpoint(2),1) ),2) )
if P.hold == 1;  hold on; end
title('core sol')
subplot(336) % outer target sol
plot(o.time, mean(o.molecule(:, max(i.grid.i_xpoint(2),1) : end ),2) )
if P.hold == 1;  hold on; end
title('outer sol')
subplot(337) % in PFR 
plot(o.time,o.extern_molecule_density(:,1));
if P.hold == 1;  hold on; end
xlabel('time (s)')
title('inner pfr')
subplot(338) % core
plot(o.time,o.core_neutral_density)
if P.hold == 1;  hold on; end
xlabel('time (s)')
title('core atom')
subplot(339) % out PFR
plot(o.time,o.extern_molecule_density(:,5));
if P.hold == 1;  hold on; end
title('outer pfr')
xlabel('time (s)')

% plasma params
figure(P.fignum + 6);

subplot(331)% DFR 1
plot(o.time,o.temperature(:,3));
ylabel('temperatures')
title('inner target')
if P.hold == 1;  hold on; end
subplot(332) % CFR
plot(o.time, mean(o.temperature(:,ceil(mean(i.grid.i_xpoint)) ),2) )
title('PLASMA omp core sol')
if P.hold == 1;  hold on; end
subplot(333) % DFR2
plot(o.time,o.temperature(:,end-2));
ylabel('temperatures')
title('outer target')
if P.hold == 1;  hold on; end
subplot(334) % inner target sol
plot(o.time,mean(o.density(:, 1: max(i.grid.i_xpoint(1),1) ),2) )
if P.hold == 1;  hold on; end
title('inner sol')
ylabel('core e density')
subplot(335) % core sol
plot(o.time, mean(o.density(:, max(i.grid.i_xpoint(1),1) : max(i.grid.i_xpoint(2),1) ),2) )
if P.hold == 1;  hold on; end
title('omp core sol')
subplot(336) % outer target sol
plot(o.time, mean(o.density(:, max(i.grid.i_xpoint(2),1) : end ),2) )
if P.hold == 1;  hold on; end
title('outer sol')
subplot(337) % in aPFR 
plot(o.time,o.Gamma_n(:,3));
if P.hold == 1;  hold on; end
xlabel('time (s)')
title('inner pfr target flux')
subplot(338) % core
plot(o.time,o.core_density)
if P.hold == 1;  hold on; end
xlabel('time (s)')
title('core plasma density')
subplot(339) % out PFR
plot(o.time,o.Gamma_n(:,end-3));
if P.hold == 1;  hold on; end
xlabel('time (s)');
title('outer pfr target flux');


%% dynamic inputs
dynstr = fields(i.dynamic);
%ncol = jet(5);
T = ceil(length(dynstr)/3);
% mylinestyles = ["-x"; "--o"; "-^";":x";"-.+"];
figure(P.fignum+ 7)
for iplot = 1:length(dynstr)
    subplot(3,T,iplot)
%     ax = gca;
%     ax.LineStyleOrder = mylinestyles;
    %ax.LineStyleCyclingMethod = "beforecolor";
    plot(o.time,i.dynamic.(dynstr{iplot})) %,'color',ncol)
    ylabel(strrep(dynstr{iplot},'_','-'));
    xlabel('t (s)'); 
    if P.hold == 1;  hold on; end
    [n,m] = size(i.dynamic.(dynstr{iplot}));
%     if m > 2
%         legend
%     end
       

end


%% plot grid
grdstr = fields(i.grid);
cntarr = 0;
for iplot = 1:length(grdstr)
 if length(i.grid.(grdstr{iplot})) > Nx-1
     if ~strcmp(grdstr{iplot},'xcb') && ~strcmp(grdstr{iplot},'x')
         cntarr = cntarr + 1;
        arrstr{cntarr} = grdstr{iplot};
     end
 end
end

T = ceil(cntarr/4);
figure(P.fignum+8);
for iplot = 1:cntarr
    subplot(4,T,iplot)
    if length(i.grid.(arrstr{iplot})) == length(x)
    plot(x,i.grid.(arrstr{iplot})) %,'color',ncol)
    elseif length(i.grid.(arrstr{iplot})) == length(xcb)
    plot(xcb,i.grid.(arrstr{iplot})) %,'color',ncol)
    end
    ylabel(strrep(arrstr{iplot},'_','-'));
    xlabel(xstr); 
    if P.hold == 1;  hold on; end
    if P.flip ==1; set(gca,'XDir','reverse');end
end

end
