function plotdiv1d_v401(output,input,varargin)
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
D.Ntime = length(output.time);
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

%% Load DIV1D Parameters
Ntime       = P.Ntime;
Nx          = length(output.X);
X           = output.X;
Xcb         = output.Xcb;
L           = max(output.X);
dxmin       = P.dxmin;
time        = output.time;
% extract quantities
density     = output.density;
temperature = output.temperature;
q_parallel  = output.q_parallel;
neutral_density = output.neutral_density;
q_X             = q_parallel(Ntime,1);
Source_v        = output.Source_v;
Gamma_mom       = output.Gamma_mom;
Gamma_n         = output.Gamma_n;
neutral_flux    = output.neutral_flux;
velocity        = output.velocity;
Source_Q        = output.Source_Q;
Source_n        = output.Source_n;
Source_neutral  = output.Source_neutral;
%B_field = input.flx_exp_prf;

% momentum
%pressure= proc.pressure; % 2.0*density.*temperature*e_charge;
%P_kin = proc.P_kin; % 2.0*density(Ntime,:).*temperature(Ntime,:)*e_charge;
%P_tot = proc.P_tot; % P_kin + Gamma_mom(Ntime,:);

P_kin = 2.0*density(Ntime,:).*temperature(Ntime,:)*e_charge;
P_tot =  P_kin + Gamma_mom(Ntime,1:end-1); 

%% General
cmap = lines;
if P.flip ==1
    X = max(X) - X;
    if isstruct(P.solps)
        solps = P.solps;
%         solps.geom.dspar = repmat(max(solps.geom.dspar),length(solps.geom.dspar(:,1)),1)-solps.geom.dspar;
    end
end  

%% Profiles
figure(P.fignum)
% temperature
subplot(2,2,3, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
plot(X(2:end-2),temperature(Ntime,2:end-2), 'LineWidth', P.LineWidth);  grid on;
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
line(1) = plot(X(2:end-2),density(Ntime,2:end-2), 'LineWidth', P.LineWidth); grid on;
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
plot(X(2:end-2),velocity(Ntime,2:end-2), 'LineWidth', P.LineWidth); grid on; 
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
plot(X(2:end-2),neutral_density(Ntime,2:end-2), 'LineWidth', P.LineWidth); grid on;
if isstruct(P.solps); hold on;
    if P.avg ==0
        plot(solps.geom.dspar(:,1),solps.data.na_correct,'--k', 'LineWidth', P.LineWidth)
    else
        plot(solps.geom.dspar(:,1),solps.data.na_correct(:,1),'--k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.na_correct(:,2),':k', 'LineWidth', P.LineWidth)
        plot(solps.geom.dspar(:,1),solps.data.na_correct(:,3),':k', 'LineWidth', P.LineWidth)
    end
end
ylabel('neutral density [1/m^3]')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end
set(gcf, 'OuterPosition',[100, 600, 500, 500]);

%% FLUXES
figure(P.fignum +1)
set(gcf, 'OuterPosition',[100, 600, 500, 500]);

% density flux
subplot(2,2,1, 'FontSize', P.FontSize); grid on;
if P.hold == 1;  hold on; end
plot(Xcb(2:end-2),Gamma_n(Ntime,2:end-2), 'LineWidth', P.LineWidth);
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
plot(Xcb(2:end-2),Gamma_mom(Ntime,2:end-2),'Color',cmap( 1,:) , 'LineWidth', P.LineWidth); hold on; grid on;
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
legend('Gamma_{mom}','P_{kin}', 'P_{tot}')
ylabel('mom flux [J/(m^3)]')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end


% energy flux
subplot(2,2,3, 'FontSize', P.FontSize)
if P.hold == 1; hold on; end
plot(Xcb(2:end-2),q_parallel(Ntime,2:end-2), 'LineWidth', P.LineWidth); grid on;
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


% neutral particle flux
subplot(2,2,4, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
plot(Xcb(2:end-2),neutral_flux(Ntime,2:end-2), 'LineWidth', P.LineWidth); grid on;
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

%% SOURCES
figure(P.fignum +2)
set(gcf, 'OuterPosition',[600, 600, 500, 500]);
subplot(2,2,1, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% density source
plot(X(2:end-2),Source_n(Ntime,2:end-2), 'LineWidth', P.LineWidth);
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
plot(X(2:end-2),Source_v(Ntime,2:end-2), 'LineWidth', P.LineWidth); grid on;
ylabel('momentum source')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end

subplot(2,2,3, 'FontSize', P.FontSize)
if P.hold == 1;  hold on; end
% plot the energy source
plot(X(2:end-2),Source_Q(Ntime,2:end-2), 'LineWidth', P.LineWidth); grid on;
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
plot(X(2:end-2),Source_neutral(Ntime,2:end-2), 'LineWidth', P.LineWidth); grid on;
ylabel('neutral source')
if P.flip ==1
    xlabel('distance to target [m]')
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream [m]')
end

        %% inputs
     figure(P.fignum+5)
     set(gcf,'OuterPosition',[600, 600, 600, 600])
       subplot(331)
       plot(output.time(2:end),input.dynamic.dyn_gas(2:end)/10^20); hold on;
       ylabel('\Gamma_{D,puff}')
        xlabel('time'); grid on; box on;
       subplot(332)
       plot(output.time(2:end),input.dynamic.dyn_nu(2:end)/10^19);hold on;
       ylabel('n_X')
       xlabel('time'); grid on; box on;
       title('inputs')
       subplot(333)
       plot(output.time(2:end),input.dynamic.dyn_qpar(2:end)/10^6);hold on;
       ylabel('qpar [MW]')
        xlabel('time'); grid on; box on;
       subplot(334)
       plot(output.time(2:end),input.dynamic.dyn_rad_los(2:end));hold on;
        xlabel('time')
        ylabel('radial los factor'); grid on; box on;
        subplot(335); 
        plot(output.time(2:end),input.dynamic.dyn_rec(2:end));hold on;
        xlabel('time');
        ylabel('recycling'); grid on; box on;
        subplot(336);
         plot(output.time(2:end),input.dynamic.dyn_red_frc(2:end));hold on;
        xlabel('time');
        ylabel('redistributed fraction'); grid on; box on;
        
        subplot(337)
        % from version 3.0.2 on there are multiple impurities as function of time
        plot(output.time(2:end),input.dynamic.dyn_imp_con(2:end,:)); hold on;
        ylabel('impurity concentration')
        xlabel('time');
        legend(num2str(input.physics.impurity_z(:)));
        grid on; box on;
        
        subplot(338)                   
        plot(output.time(2:end),input.dynamic.dyn_nb(2:end));hold on;
        xlabel('time');
        ylabel('neutral background'); grid on; box on;
        
         
        subplot(339)
        plot(input.physics.l-output.X,input.grid.b_field);hold on;
        plot(input.physics.l-output.X,input.grid.gas_puff_profile);
        plot(input.physics.l-output.X,input.grid.core_source_profile);
        ylabel('a.u.');
        xlabel('distance to target');
        legend('B_field','puff profile','core profile');
        set(gca,'XDir','reverse'); grid on; box on;

      
       
  
end
