function plot_div1d_output(div1doutput,div1dinput,proc,varargin)
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
D.Ntime = length(div1doutput.time);
D.flip = 1;
D.FontSize = 11;%2;
D.fignum = 16;
D.LineWidth = 1;
D.hold = 0;
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

%% Load DIV1D Parameters
Ntime       = P.Ntime;
Nx          = length(div1doutput.X);
X           = div1doutput.X;
Xcb         = div1doutput.Xcb;
L           = max(div1doutput.X);
dxmin       = P.dxmin;
time = div1doutput.time;
% extract quantities
density     = div1doutput.density;
temperature = div1doutput.temperature;
q_parallel  = div1doutput.q_parallel;
neutral_density = div1doutput.neutral_density;
q_X         = q_parallel(Ntime,1);
Source_v    = div1doutput.Source_v;
Gamma_mom   = div1doutput.Gamma_mom;
Gamma_n     = div1doutput.Gamma_n;
neutral_flux = div1doutput.neutral_flux;
velocity    = div1doutput.velocity;
Source_Q    = div1doutput.Source_Q;
Source_n    = div1doutput.Source_n;
Source_neutral = div1doutput.Source_neutral;
Source_v = div1doutput.Source_v;
B_field = div1dinput.grid.b_field.*div1dinput.grid.b_trans; %div1dinput.flx_exp_prf;

% momentum
pressure= proc.pressure; % 2.0*density.*temperature*e_charge;
P_kin = proc.P_kin; % 2.0*density(Ntime,:).*temperature(Ntime,:)*e_charge;
P_tot = proc.P_tot; % P_kin + Gamma_mom(Ntime,:); %.*B_field;


%% General
cmap = lines;
if P.flip ==1
    X = max(Xcb) - X;
    Xcb = max(Xcb)- Xcb;
    if isstruct(P.solps)
        solps = P.solps;
%         solps.geom.dspar = repmat(max(solps.geom.dspar),length(solps.geom.dspar(:,1)),1)-solps.geom.dspar;
    end
end  

%% equation terms
figure(P.fignum -1)
subplot(2,2,3) % energy
plot(X, Source_Q(Ntime,:)); hold on;
 plot(X, proc.dqdt);
%plot(X, proc.ene_bdqdx+ Source_Q(Ntime,:));
plot(X, proc.ene_bdqdx,'--'); 
plot(X, proc.ene_vdpdx,'--'); hold off;
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
ylabel('energy terms');
legend('Q','dqdt','bdqdx','vdpdx');
title('check equations')

subplot(2,2,1) % density
plot(X, Source_n(Ntime,:)); hold on;
plot(X, proc.dndt);
plot(X, proc.den_rhs,'--'); hold off;
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
ylabel('density terms');
legend('S n','dndt','den rhs');

subplot(2,2,2) % momentum
plot(X, Source_v(Ntime,:)); hold on;
plot(X, proc.dmomdt); 
plot(X, proc.mom_bdgamdx,'--'); 
plot(X, proc.mom_dpdx,'--'); 
hold off;
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
ylabel('mom terms');
legend('S v','dvdt','bdgamdx','dpdx');

subplot(2,2,4) % neutral
plot(X, Source_neutral(Ntime,:)); hold on;
plot(X, proc.dnndt);
plot(X, proc.neu_diffdnndx,'--'); hold off;
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
ylabel('neutral terms');
legend('S neut','dnndt','diffdnndx');

%% Profiles
figure(P.fignum)
% temperature
subplot(2,2,3, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
plot(X,temperature(Ntime,:), 'LineWidth', P.LineWidth); hold on;  grid on;
if isstruct(P.D2PM)
if P.flip ==1
    line(1) = plot(flip([0,D2PM.L]),[D2PM.T_X(Ntime,1) D2PM.T_L(Ntime,1)],'+','Color',cmap(1,:),'LineWidth', P.LineWidth);
    line(2) = plot(flip([0,D2PM.L]),[D2PM.T_X(Ntime,2) D2PM.T_L(Ntime,2)],'x','Color',cmap(1,:),'LineWidth', P.LineWidth);
else
    line(1) = plot([0,D2PM.L],[D2PM.T_X(Ntime,1) D2PM.T_L(Ntime,1)],'+','Color',cmap(1,:),'LineWidth', P.LineWidth);
    line(2) = plot([0,D2PM.L],[D2PM.T_X(Ntime,2) D2PM.T_L(Ntime,2)],'x','Color',cmap(1,:),'LineWidth', P.LineWidth);
end
end
if isstruct(P.solps);   hold on;
    ticol = [0.6350, 0.0780, 0.1840];
    if P.avg ==0
        line(3) = plot(solps.geom.dspar,solps.data.te,'--k', 'LineWidth', P.LineWidth);
        line(4) = plot(solps.geom.dspar,solps.data.ti,'--','Color',ticol, 'LineWidth', P.LineWidth);
    else % also plot the range that was averaged
        line(3) = plot(solps.geom.dspar(:,1),solps.data.te(:,1),'--k', 'LineWidth', P.LineWidth);   % mean
        plot(solps.geom.dspar(:,1),solps.data.te(:,2),':k', 'LineWidth', P.LineWidth) % minimum
        plot(solps.geom.dspar(:,1),solps.data.te(:,3),':k', 'LineWidth', P.LineWidth) % maximum
        line(4) =plot(solps.geom.dspar(:,1),solps.data.ti(:,1),'--','Color',ticol, 'LineWidth', P.LineWidth);   % mean ;
        plot(solps.geom.dspar(:,1),solps.data.ti(:,2),':', 'Color',ticol,'LineWidth', P.LineWidth) % minimum
        plot(solps.geom.dspar(:,1),solps.data.ti(:,3),':','Color',ticol, 'LineWidth', P.LineWidth) % maximum
    end
    if isstruct(P.S2PM)
    if P.flip ==1
        plot(flip([0,S2PM.L]),[S2PM.T_X(1) S2PM.T_L(1)],'+k','LineWidth', P.LineWidth);
        plot(flip([0,S2PM.L]),[S2PM.T_X(2) S2PM.T_L(2)],'xk','LineWidth', P.LineWidth);
    else
        plot([0,S2PM.L],[S2PM.T_X(1) S2PM.T_L(1)],'+k','LineWidth', P.LineWidth);
        plot([0,S2PM.L],[S2PM.T_X(2) S2PM.T_L(2)],'xk','LineWidth', P.LineWidth);
    end
    legend(line,{'2PM','K-R','T_e','T_i'});
    end
    legend(line(3:4),{'T_e','T_i'});
end
ylabel('temperature ([V]');

if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
title(P.Dtitle);

% density
subplot(2,2,1, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
line(1) = plot(X,density(Ntime,:), 'LineWidth', P.LineWidth); hold on; grid on;
if isstruct(P.D2PM)
if P.flip ==1
    plot(0,D2PM.n_L(Ntime,1),'+','Color',cmap(1,:),'LineWidth', P.LineWidth);
    plot(0,D2PM.n_L(Ntime,2),'x','Color',cmap(1,:),'LineWidth', P.LineWidth);
else
    plot(D2PM.L,D2PM.n_L(Ntime,1),'+','Color',cmap(1,:),'LineWidth', P.LineWidth);
    plot(D2PM.L,D2PM.n_L(Ntime,2),'x','Color',cmap(1,:),'LineWidth', P.LineWidth);
end
end
if isstruct(P.solps);   hold on;
    if P.avg ==0
        line(2) = plot(solps.geom.dspar,solps.data.ne,'--k', 'LineWidth', P.LineWidth);
    else
        line(2) = plot(solps.geom.dspar(:,1),solps.data.ne(:,1),'--k', 'LineWidth', P.LineWidth);
        plot(solps.geom.dspar(:,1),solps.data.ne(:,2),':k', 'LineWidth', P.LineWidth);
        plot(solps.geom.dspar(:,1),solps.data.ne(:,3),':k', 'LineWidth', P.LineWidth);
    end
    if isstruct(P.S2PM)
    if P.flip ==1
        plot(0,S2PM.n_L(1),'+k','LineWidth', P.LineWidth);
        plot(0,S2PM.n_L(2),'xk','LineWidth', P.LineWidth);
    else
        plot(S2PM.L,S2PM.n_L(1),'+k','LineWidth', P.LineWidth);
        plot(S2PM.L,S2PM.n_L(2),'xk','LineWidth', P.LineWidth);
    end
    end
end
legend(line,{'DIV1D','SOLPS'},'location','best');
ylabel('density [/m^3]');
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
if isstruct(P.solps)
    ylim([0 max(solps.data.ne(:,3))*1.05]);
end
title(P.title);

% velocity
subplot(2,2,2, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
plot(X,velocity(Ntime,:), 'LineWidth', P.LineWidth); grid on; hold on;
if isstruct(P.solps);   hold on;
    % species dependent(D0, D+1, C0, C+1:6)
    if P.avg == 0
        plot(solps.geom.dspar,-1*solps.datat.ua(:,2),'--k', 'LineWidth', P.LineWidth)
    else
        plot(solps.geom.dspar(:,1),solps.data.ua(:,1,2),'--k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.ua(:,2,2),':k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.ua(:,3,2),':k', 'LineWidth', P.LineWidth)
    end
end
ylabel('parallel velocity [m/s]')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end

% neutral density
subplot(2,2,4, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
plot(X,neutral_density(Ntime,:), 'LineWidth', P.LineWidth); grid on;
if isstruct(P.solps); hold on;
    if P.avg ==0
        plot(solps.geom.dspar(:,1),solps.data.na_correct,'--k', 'LineWidth', P.LineWidth)
    else
        plot(solps.geom.dspar(:,1),solps.data.na_correct(:,1),'--k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.na_correct(:,2),':k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.na_correct(:,3),':k', 'LineWidth', P.LineWidth)
    end
end
ylabel('neutral density [1/m3]','FontSize', P.FontSize)
if P.flip ==1
    xlabel('distance to target [m]','FontSize', P.FontSize)
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]','FontSize', P.FontSize)
end


set(gcf, 'OuterPosition',[100, 600, 500, 500]);

%% FLUXES
figure(P.fignum +1)
set(gcf, 'OuterPosition',[100, 100, 500, 500]);

% density flux
subplot(2,2,1, 'FontSize', P.FontSize); grid on;
if P.hold == 1;  hold on; end
plot(Xcb,Gamma_n(Ntime,:), 'LineWidth', P.LineWidth);
if isstruct(P.solps);   hold on;
    %   species dependent(D0, D+1, C0, C+1:6)
    if P.avg == 0
        plot(solps.geom.dspar(1),solps.data.ne.*solps.data.ua(:,2),'--k', 'LineWidth', P.LineWidth)
    else
        plot(solps.geom.dspar(:,1),solps.data.ne(:,1).*solps.data.ua(:,1,2),'--k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.ne(:,2).*solps.data.ua(:,2,2),':k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.ne(:,3).*solps.data.ua(:,3,2),':k', 'LineWidth', P.LineWidth)
    end
end
ylabel('part flux [1/(m^2s)')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;

% momentum flux
subplot(2,2,2, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
plot(Xcb,Gamma_mom(Ntime,:),'Color',cmap( 1,:) , 'LineWidth', P.LineWidth); hold on;
plot(X,P_kin,'Color',cmap(3+2,:) , 'LineWidth', P.LineWidth);
plot(X,P_tot,'Color',cmap(3+3,:) , 'LineWidth', P.LineWidth);
if isstruct(P.solps); hold on;
    if P.avg ==0
        plot(solps.geom.dspar,solps.data.pdyn,'--','Color',cmap(1,:) , 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar,solps.data.pstat,'--','Color',cmap(2,:), 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar,solps.data.ptot,'--','Color',cmap(3,:), 'LineWidth', P.LineWidth)
    else
        plot(solps.geom.dspar(:,1),solps.data.pdyn(:,1),'--','Color',cmap(1,:) , 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.pstat(:,1),'--','Color',cmap(2,:), 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.ptot(:,1),'--','Color',cmap(3,:), 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.pdyn(:,2),':','Color',cmap(1,:) , 'LineWidth', P.LineWidth-.5)
        plot(solps.geom.dspar(:,1),solps.data.pstat(:,2),':','Color',cmap(2,:), 'LineWidth', P.LineWidth-.5)
        plot(solps.geom.dspar(:,1),solps.data.ptot(:,2),':','Color',cmap(3,:), 'LineWidth', P.LineWidth-.5)
        plot(solps.geom.dspar(:,1),solps.data.pdyn(:,3),':','Color',cmap(1,:) , 'LineWidth', P.LineWidth-.5)
        plot(solps.geom.dspar(:,1),solps.data.pstat(:,3),':','Color',cmap(2,:), 'LineWidth', P.LineWidth-.5)
        plot(solps.geom.dspar(:,1),solps.data.ptot(:,3),':','Color',cmap(3,:), 'LineWidth', P.LineWidth-.5)
    end
end
legend('Gamma mom','P kin', 'P tot')
ylabel('mom flux [J/(m^3)]')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;

% energy flux
subplot(2,2,3, 'FontSize', P.FontSize)
if P.hold == 1; hold on; end
plot(Xcb,q_parallel(Ntime,:), 'LineWidth', P.LineWidth);
if isstruct(P.solps);   hold on;
    if P.avg ==0
        plot(solps.geom.dspar,abs(solps.data.qpar),'--k', 'LineWidth', P.LineWidth)
    else
        plot(solps.geom.dspar(:,1),abs(solps.data.qpar(:,1)),'--k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),abs(solps.data.qpar(:,2)),':k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),abs(solps.data.qpar(:,3)),':k', 'LineWidth', P.LineWidth)
    end
end
ylabel('energy flux [W/m^2]')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;

% neutral particle flux
subplot(2,2,4, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
plot(Xcb,neutral_flux(Ntime,:), 'LineWidth', P.LineWidth);
if isstruct(P.solps);   hold on;
    %   species dependent(D0, D+1, C0, C+1:6)
    if P.avg == 0
        plot(solps.geom.dspar,-1*solps.data.na_correct(:,1).*solps.data.ua(:,1),'--k', 'LineWidth', P.LineWidth)
    else
        plot(solps.geom.dspar,-1*solps.data.na_correct(:,1,1).*solps.data.ua(:,1,1),'--k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar,-1*solps.data.na_correct(:,2,1).*solps.data.ua(:,2,1),':k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar,-1*solps.data.na_correct(:,2,1).*solps.data.ua(:,2,1),':k', 'LineWidth', P.LineWidth)
    end
end
ylabel('neutral flux [1/(m^2s)]')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;

%% SOURCES
figure(P.fignum +3)
set(gcf, 'OuterPosition',[100, 100, 500, 500]);
subplot(2,2,1, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% density source
plot(X,Source_n(Ntime,:), 'LineWidth', P.LineWidth);
title('sources')
ylabel('density source')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;
subplot(2,2,2, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% plot the momentum source
plot(X,Source_v(Ntime,:), 'LineWidth', P.LineWidth);
ylabel('momentum source')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;
subplot(2,2,3, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% plot the energy source
plot(X(1:end-3),Source_Q(Ntime,1:end-3), 'LineWidth', P.LineWidth); hold on;
% plot(X(1:end-3),Source_Qtot(1:end-3),'LineWidth', P.LineWidth);
if isstruct(P.solps)
    hold on;
    clear line
    if P.avg ==0
        line(1)=  plot(solps.geom.dspar,-solps.data.totrad_c,'--','Color',cmap(3,:) , 'LineWidth', P.LineWidth);
    else
        line(1)=  plot(solps.geom.dspar(:,1),-sum(solps.data.totrad_h(:,1,:),3),'--','Color',cmap(3,:) , 'LineWidth', P.LineWidth);
        line(2)=  plot(solps.geom.dspar(:,1),-sum(solps.data.totrad_c(:,1,:),3),'--','Color',cmap(4,:) , 'LineWidth', P.LineWidth);
        line(3)=  plot(solps.geom.dspar(:,1),-sum(solps.data.totrad2(:,1,:),3),'--','Color',cmap(5,:) , 'LineWidth', P.LineWidth);
        plot(solps.geom.dspar(:,1),-sum(solps.data.totrad_h(:,2,:),3),':','Color',cmap(3,:) , 'LineWidth', P.LineWidth);
        plot(solps.geom.dspar(:,1),-sum(solps.data.totrad_c(:,2,:),3),':','Color',cmap(4,:) , 'LineWidth', P.LineWidth);
        plot(solps.geom.dspar(:,1),-sum(solps.data.totrad2(:,2,:),3),':','Color',cmap(5,:) , 'LineWidth', P.LineWidth);
        plot(solps.geom.dspar(:,1),-sum(solps.data.totrad_h(:,3,:),3),':','Color',cmap(3,:) , 'LineWidth', P.LineWidth);
        plot(solps.geom.dspar(:,1),-sum(solps.data.totrad_c(:,3,:),3),':','Color',cmap(4,:) , 'LineWidth', P.LineWidth);
        plot(solps.geom.dspar(:,1),-sum(solps.data.totrad2(:,3,:),3),':','Color',cmap(5,:) , 'LineWidth', P.LineWidth);
    end
    legend(line(1:3),{'rad H','rad C','rad tot'},'location','best');
end

ylabel('energy source [SOLPS = W/m3]')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
% ylim([min(Source_Qtot(1:end-3)) 0]);
grid on;
subplot(2,2,4, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% plot the neutral source
plot(X,Source_neutral(Ntime,:), 'LineWidth', P.LineWidth);
ylabel('neutral source')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
grid on;

%% detailed SOURCES 
% % TODO plot the detailed sources against the outputs of DIV1D.
% figure(P.fignum +7)
% set(gcf, 'OuterPosition',[600, 600, 500, 500]);
% title('detailed sources')
% % density sources
% subplot(2,2,1, 'FontSize', P.FontSize)
% plot(X,Source_n(Ntime,:), 'LineWidth', P.LineWidth); hold on;
% hold off;
% ylabel('density source')
% if P.flip ==1
%     xlabel('distance to target [m]')
%     set(gca,'XDir','reverse');
% else
%     xlabel('distance from upstream [m]')
% end
% grid on;
% 
% % momentum sources
% subplot(2,2,2, 'FontSize', P.FontSize)
% plot(X,Source_v(Ntime,:), 'LineWidth', P.LineWidth); hold on;
% hold off
% ylabel('momentum source')
% if P.flip ==1
%     xlabel('distance to target [m]')
%     set(gca,'XDir','reverse');
% else
%     xlabel('distance from upstream [m]')
% end
% grid on;
% 
% %  energy sources
% subplot(2,2,3, 'FontSize', P.FontSize)
% plot(X(1:end-3),Source_Q(Ntime,1:end-3), 'LineWidth', P.LineWidth); hold on;
% hold off;
% ylabel('energy sources')
% if P.flip ==1
%     xlabel('distance to target [m]')
%     set(gca,'XDir','reverse');
% else
%     xlabel('distance from upstream [m]')
% end
% grid on;
% 
% % neutral sources
% subplot(2,2,4, 'FontSize', P.FontSize)
% plot(X,Source_neutral(Ntime,:), 'LineWidth', P.LineWidth); hold on;
% hold off;
% ylabel('neutral source')
% if P.flip ==1
%     xlabel('distance to target [m]')
%     set(gca,'XDir','reverse');
% else
%     xlabel('distance from upstream [m]')
% end
% grid on;
% 

        %% inputs
     figure(P.fignum+5)
     set(gcf,'OuterPosition',[1100, 200, 600, 600])
       subplot(331)
       plot(div1doutput.time(2:end),div1dinput.dynamic.dyn_gas(2:end)/10^20); hold on;
       ylabel('dyn gas')
        xlabel('time'); grid on; box on;
       subplot(332)
       plot(div1doutput.time(2:end),div1dinput.dynamic.dyn_nu(2:end)/10^19);hold on;
       ylabel('dyn n X')
       xlabel('time'); grid on; box on;
       title('inputs')
       subplot(333)
       plot(div1doutput.time(2:end),div1dinput.dynamic.dyn_qpar(2:end)/10^6);hold on;
       ylabel('qpar [MW]')
        xlabel('time'); grid on; box on;
       subplot(334)
%        plot(div1doutput.time(2:end),div1dinput.dynamic.dyn_rad_los(2:end));hold on;
        %xlabel('time')
       % ylabel('radial los factor'); grid on; box on;
        subplot(335); 
        plot(div1doutput.time(2:end),div1dinput.dynamic.dyn_rec(2:end));hold on;
        xlabel('time');
        ylabel('recycling'); grid on; box on;
        subplot(336);
         plot(div1doutput.time(2:end),div1dinput.dynamic.dyn_red_frc(2:end));hold on;
        xlabel('time');
        ylabel('redistributed fraction'); grid on; box on;
        
        subplot(337)
%         try
%         plot(div1dinput.L-div1doutput.X,div1dinput.grid.car_con_prf);hold on;
%         ylabel('carbon profile');
%         xlabel('distance to target');
%         set(gca,'XDir','reverse'); 
%         catch % from version 3.0.2 on there are multiple impurities as function of time
         
            plot(div1doutput.time(2:end),div1dinput.dynamic.dyn_imp_con(2:end,:)); hold on;
            ylabel('impurity concentration')
            xlabel('time');
            legend(num2str(div1dinput.physics.impurity_z(:)));
%         end
        grid on; box on;
          subplot(338)
         
          
        % this might become multiple in the latest version
        plot(div1doutput.time(2:end),div1dinput.dynamic.dyn_nb(2:end));hold on; 
        xlabel('time');
        ylabel('neutral background'); grid on; box on;
                 
%         plot(div1dinput.L-div1doutput.X,div1dinput.gas_puf_prf);hold on;
%         ylabel('gas puff profile');
%         xlabel('distance to target');
%         set(gca,'XDir','reverse'); grid on; box on;
        subplot(339)
        try
            plot(div1dinput.physics.l-div1doutput.X,div1dinput.grid.b_field);hold on;
        catch
            plot(div1dinput.physics.l-div1doutput.X,div1dinput.grid.flx_exp_prf);hold on; 
        end
        plot(div1dinput.physics.l-div1doutput.X,div1dinput.grid.gas_puff_profile);
        ylabel('a.u.');
        xlabel('distance to target');
        legend('B field','gas puff profile');
        set(gca,'XDir','reverse'); grid on; box on;

%         figure
       
  
end
