function [h] = plot_div1d_profiles_v600(o, i, varargin)
% function used to plot the div1d profiles (together with other profiles)

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

D.gamma = 6;
D.dxmin = 0.05; 
D.xlimits = [0  17.8];
D.Ntime = length(o.time);
D.solps = 0;
D.input = 0; % div1d input struct
D.version = 'v2.0.0';


% -- plotting ---
D.plotoutput = 0;
D.ploterror = 0;
D.plotinput = 0;
D.timing = 0;
D.generalposition = [100 800 200 250];
D.print = 0;
D.printerror =1;
D.plot = 0;
D.figtight = 0;
D.plot2PM = 0;
D.allpossiblefigures = 0;
D.seperatefigs =0;
D.colorindex = 0;
D.fignum = 16;
D.LineWidth = 1;
D.FontSize = 11;
D.flip = 1;
D.title = '' ;
D.Dtitle = '' ;
D.save = 0;
D.cratio = 0; 
D.hold = 0;
D.fR = 0;
D.R = 1.0;
D.avg = 1;
D.FE = 1.0;
D.tn = 1;
D.nb    = [0 0 0 0 0];
D.xticks = [0 1 2 3 4 5 6 7 8 9 10 15];
D.xp = 15.6;
D.ytic_q = 0; % heat flux
D.ylim_q = 0;
D.ytic_t = 0; % temperature
D.ylim_t = 0;
D.ytic_n = 0; % electron density
D.ylim_n = 0;
D.ytic_v = 0; % velocity
D.ylim_v = 0;
D.ytic_m = 0; % molecules
D.ylim_m = 0;
D.ytic_nv = 0;
D.ylim_nv = 0;
D.ytic_a = 0;
D.ylim_a = 0;
D.ytic_f = 0;
D.ylim_f = 0; 
D.tf_form = [0.2 0.04 0.06 0.04]; % dx_left dx_right, dy_bot, dy_top
D.qlegloc = 'best';
D.tlegloc = 'best'; 
D.implabel = '';
D.plotimpurities = 1;
D.plotbackground = 1;
D.pxline = 1;
D.i_xpoint = [0 0];
P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end
P.tf_form = min(1,max(P.tf_form,0)); % defined inside the figure.
plimp = P.plotimpurities; plimp = 0; % plimp -1;
Ntime       = P.Ntime;
Nx          = length(o.X);
X           = o.X;
Xcb         = o.Xcb;
L           = max(o.X);
dxmin       = P.dxmin;

B_field = i.grid.b_field';

if isstruct(P.solps)
solps = P.solps;
end

if P.flip ==1
    X = max(X) - X;
    Xcb = max(Xcb) - Xcb;
    if isstruct(P.solps)
        solps = P.solps;
%          solps.geom.dspar = repmat(max(solps.geom.dspar),length(solps.geom.dspar(:,1)),1)-solps.geom.dspar;
    end
end  

    % reload figure
     %figure(P.fignum+4);
    fig = gcf;
    keeplabel = true; 
    if fig.Number == P.fignum+4 
        isfig = true; % there is already a figure, no need to redo all

    else
        isfig = false;
        figure(P.fignum+4);
    end
    generalposition  = P.generalposition;
    
    % energy flux
    
    if P.figtight ==0
        position = generalposition + [400 600 0 0 ];
    else
        position = generalposition; %[100 200 300 600];% [x y width height]        
        ax = struct; 
        x1 = P.tf_form(1); % .16;%0.25;
        x2 = P.tf_form(2); %x1;
        %xp = 5.58;

        y1 = P.tf_form(3); % 0.06  %0.96-yo;
        y2 = P.tf_form(4); % y1; 
        w = 1-x2-x1;
        h = 1-y2-y1;
        dy = h/7;
        qpos    = [x1, y1+dy*6,w,dy ];
        tpos    = [x1, y1+dy*5,w,dy ];
        nepos   = [x1, y1+dy*4,w,dy ];
        upos    = [x1, y1+dy*3,w,dy ];
        ndpos   = [x1, y1+dy*2,w,dy ];
        nvpos   = [x1, y1+dy*1,w,dy ];
        mpos    = [x1, y1     ,w,dy ];
    end
    if ~isfig ; set(gcf, 'Position',position); end
    if P.figtight ==1
        if ~isfig 
            ax.q = axes('Position',qpos, 'FontSize', P.FontSize-1); 
        else
            ax.q = fig.Children(7+plimp);  
            fig.CurrentAxes = ax.q; 
        end 
       
        box on; 
    end

    
    if P.hold == 1;  hold on; end
   
    qparline(1) =  plot(Xcb,o.q_parallel(Ntime,:)/10^6, 'LineWidth', P.LineWidth);
 if P.timing == 1
        tmpsring = strcat('t=',num2str(round(o.time(Ntime)*1000,1,"decimals")),' [ms]');
        text(15,10,tmpsring);
    end


    if P.pxline ==1
    if i.grid.i_xpoint(1) > 1;  xline(X(i.grid.i_xpoint(1)),'-'); end
    if i.grid.i_xpoint(2) > 1;  xline(X(i.grid.i_xpoint(2)),'-'); end
    end
    if isstruct(P.solps);   hold on;
        if P.avg ==0
            plot(solps.geom.dspar,abs(solps.data.qpar)/10^6,'--k', 'LineWidth', P.LineWidth)
        else
            tmpX = [solps.geom.dspar(1,1); solps.geom.dspar(:,1); flip(solps.geom.dspar(:,1))];
            tmpY = [solps.data.qpar(1,1); solps.data.qpar(:,2);  flip(solps.data.qpar(:,3))]/10^6;
            patch(tmpX,tmpY,[0 0 0],'FaceAlpha',0.3,'LineStyle','none')
            qparline(2) =  plot(solps.geom.dspar(:,1),abs(solps.data.qpar(:,1))/10^6,'--k', 'LineWidth', P.LineWidth);
            % the errorbar
        end
    end

    if  keeplabel; qlab = ylabel('$q_{\|}$ (MW m$^{-2}$)','FontSize',P.FontSize,'interpreter','latex'); end 

    if isstruct(P.solps)
        legend(qparline,{'DIV1D','SOLPS'},'FontSize',P.FontSize-2,'interpreter','latex','location',P.qlegloc);
    end

    if P.flip ==1
        if P.figtight ==0;    xlabel('distance to target (m)','FontSize',P.FontSize,'interpreter','latex');
        else
            xticklabels({''}); xticks(P.xticks);
        end
        set(gca,'XDir','reverse');
    else
        xlabel('distance from upstream (m)','FontSize',P.FontSize,'interpreter','latex')
    end
    grid on;
    
    if ~isempty(P.title); if ~isfig; title(P.title,'FontSize',P.FontSize,'interpreter','latex'); end; end
    xlim(P.xlimits);
    if max(P.ylim_q) ~=0;    ylim(P.ylim_q); end
    if max(P.ytic_q) ~=0;  yticks(P.ytic_q); end
    
    % temperature
    if P.figtight ==0
        figure(P.fignum+5);
        position = generalposition + [400 300 0 0 ];
        set(gcf, 'OuterPosition',position);
    else
        if ~isfig; ax.t = axes('Position',tpos, 'FontSize', P.FontSize-1); end %else; ax.t = fig.Children(6+plimp);  fig.CurrentAxes = ax.t; end; box on;
    end
    if P.hold == 1;  hold on; end
    plot(X,o.temperature(Ntime,:), 'LineWidth', P.LineWidth);  grid on;
    if P.pxline ==1
    if i.grid.i_xpoint(1) > 1;  xline(X(i.grid.i_xpoint(1)),'-'); end
    if i.grid.i_xpoint(2) > 1;  xline(X(i.grid.i_xpoint(2)),'-'); end
    end
    if P.plot2PM ==1
        if P.flip ==1
            line(1) = plot(flip([0,D2PMwivi.L]),[D2PMwivi.T_X(Ntime,1) D2PMwivi.T_L(Ntime,1)],'+','Color',cmap(1,:),'LineWidth', P.LineWidth);
            line(2) = plot(flip([0,D2PMwivi.L]),[D2PMwivi.T_X(Ntime,2) D2PMwivi.T_L(Ntime,2)],'x','Color',cmap(1,:),'LineWidth', P.LineWidth);
        else
            line(1) = plot([0,D2PMwivi.L],[D2PMwivi.T_X(Ntime,1) D2PMwivi.T_L(Ntime,1)],'+','Color',cmap(1,:),'LineWidth', P.LineWidth);
            line(2) = plot([0,D2PMwivi.L],[D2PMwivi.T_X(Ntime,2) D2PMwivi.T_L(Ntime,2)],'x','Color',cmap(1,:),'LineWidth', P.LineWidth);
        end
    end
    if isstruct(P.solps);   hold on;
        ticol = [0.6350, 0.0780, 0.1840];
        if P.avg ==0
            line(3) = plot(solps.geom.dspar,solps.data.te,'--k', 'LineWidth', P.LineWidth);
            line(4) = plot(solps.geom.dspar,solps.data.ti,'--','Color',ticol, 'LineWidth', P.LineWidth);
        else % also plot the range that was averaged
            line(3) = plot(solps.geom.dspar(:,1),solps.data.te(:,1),'--k', 'LineWidth', P.LineWidth);   % mean
            tmpX = [solps.geom.dspar(1,1); solps.geom.dspar(:,1); flip(solps.geom.dspar(:,1))];
            tmpY = [solps.data.te(1,1); solps.data.te(:,2);  flip(solps.data.te(:,3))];
            patch(tmpX,tmpY,[0 0 0],'FaceAlpha',0.3,'LineStyle','none') % the errorbar
            line(4) =plot(solps.geom.dspar(:,1),solps.data.ti(:,1),'--','Color',ticol, 'LineWidth', P.LineWidth);   % mean ;
            tmpX = [solps.geom.dspar(1,1); solps.geom.dspar(:,1); flip(solps.geom.dspar(:,1))];
            tmpY = [solps.data.ti(1,1); solps.data.ti(:,2);  flip(solps.data.ti(:,3))];
            patch(tmpX,tmpY,ticol,'FaceAlpha',0.3,'LineStyle','none') % the errorbar
        end
        if P.plot2PM==1
            if P.flip ==1
                plot(flip([0,S2PM.L]),[S2PM.T_X(1) S2PM.T_L(1)],'+k','LineWidth', P.LineWidth);
                plot(flip([0,S2PM.L]),[S2PM.T_X(2) S2PM.T_L(2)],'xk','LineWidth', P.LineWidth);
            else
                plot([0,S2PM.L],[S2PM.T_X(1) S2PM.T_L(1)],'+k','LineWidth', P.LineWidth);
                plot([0,S2PM.L],[S2PM.T_X(2) S2PM.T_L(2)],'xk','LineWidth', P.LineWidth);
            end
        end
    end
    if P.plot2PM ==1
        legend(line,{'2PM','K-R','$T_{\mathrm{e}}$','$T_{\mathrm{i}}$'},'FontSize',P.FontSize-2,'interpreter','latex'); %,'interpreter','latex'
    else
        if isstruct(P.solps)
            legend(line(3:4),{'$T_{\mathrm{e}}$','$T_{\mathrm{i}}$'},'FontSize',P.FontSize-2,'interpreter','latex','location',P.tlegloc);   %,'interpreter','latex'
        end
    end
    if  keeplabel; tlab = ylabel('$T$ (eV)','FontSize',P.FontSize,'interpreter','latex'); end
     if max(P.ylim_t) ~=0;    ylim(P.ylim_t); end
     if max(P.ytic_t) ~=0;  yticks(P.ytic_t); end
    

    if P.flip ==1
        if P.figtight ==0;    xlabel('distance to target (m)','FontSize',P.FontSize,'interpreter','latex'); %,'interpreter','latex',
        else
             xticklabels({''}); xticks(P.xticks);
        end
        set(gca,'XDir','reverse');
    else
        xlabel('distance from upstream (m)','FontSize',P.FontSize,'interpreter','latex')
    end
    xlim(P.xlimits);


    % density electron
    if P.figtight ==0
        figure(P.fignum+6);
        position = generalposition + [400 0 0 0 ];
        set(gcf, 'OuterPosition',position);
    else
       if ~isfig; ax.ne = axes('Position',nepos, 'FontSize', P.FontSize-1); else; ax.ne = fig.Children(5+plimp);  fig.CurrentAxes = ax.ne; end; ax.ne; box on;
    end
    if P.hold == 1;  hold on; end
    plot(X,o.density(Ntime,:), 'LineWidth', P.LineWidth); grid on; 
    if P.pxline == 1
    if i.grid.i_xpoint(1) > 1;  xline(X(i.grid.i_xpoint(1)),'-'); end
    if i.grid.i_xpoint(2) > 1;  xline(X(i.grid.i_xpoint(2)),'-'); end
    end
    if P.plot2PM ==1
        if P.flip ==1
            plot(0,D2PMwivi.n_L(Ntime,1),'+','Color',cmap(1,:),'LineWidth', P.LineWidth);
            plot(0,D2PMwivi.n_L(Ntime,2),'x','Color',cmap(1,:),'LineWidth', P.LineWidth);
        else
            plot(D2PMwivi.L,D2PMwivi.n_L(Ntime,1),'+','Color',cmap(1,:),'LineWidth', P.LineWidth);
            plot(D2PMwivi.L,D2PMwivi.n_L(Ntime,2),'x','Color',cmap(1,:),'LineWidth', P.LineWidth);
        end
    end
    if isstruct(P.solps);   hold on;
        if P.avg ==0
            plot(solps.geom.dspar(:,1),solps.data.ne,'--k', 'LineWidth', P.LineWidth)
        else
            plot(solps.geom.dspar(:,1),solps.data.ne(:,1),'--k', 'LineWidth', P.LineWidth)
            tmpX = [solps.geom.dspar(1,1); solps.geom.dspar(:,1); flip(solps.geom.dspar(:,1))];
            tmpY = [solps.data.ne(1,1); solps.data.ne(:,2);  flip(solps.data.ne(:,3))];
            patch(tmpX,tmpY,[0 0 0],'FaceAlpha',0.3,'LineStyle','none') % the errorbar
        end
        if P.plot2PM ==1
            if P.flip ==1
                plot(0,S2PM.n_L(1),'+k','LineWidth', P.LineWidth);
                plot(0,S2PM.n_L(2),'xk','LineWidth', P.LineWidth);
            else
                plot(S2PM.L,S2PM.n_L(1),'+k','LineWidth', P.LineWidth);
                plot(S2PM.L,S2PM.n_L(2),'xk','LineWidth', P.LineWidth);
            end
        end
    end
    if  keeplabel; nlab = ylabel('$n_{\mathrm{e}}$ (m$^{-3}$)','FontSize',P.FontSize,'interpreter','latex'); end 
    set(gca,'YScale','log')
    if P.ylim_n ~=0;    ylim(P.ylim_n); end
    if P.ytic_n ~=0;  yticks(P.ytic_n); end
    if P.flip ==1

    if P.figtight ==0;    xlabel('distance to target (m)','FontSize',P.FontSize,'interpreter','latex'); 
        else
             xticklabels({''}); xticks(P.xticks);
        end
        set(gca,'XDir','reverse');
    else
        xlabel('distance from upstream (m)','interpreter','latex','FontSize',P.FontSize)
    end
    xlim(P.xlimits);

    % velocity
    if P.figtight ==0
        figure(P.fignum+7);
        position = generalposition + [0 600 0 0 ];
        set(gcf, 'OuterPosition',position);
    else
       if ~isfig; ax.u = axes('Position',upos, 'FontSize', P.FontSize-1); else; ax.u = fig.Children(4+plimp);  fig.CurrentAxes = ax.u; end; ax.u; box on;
    end
    if P.hold == 1;  hold on; end
    plot(X,o.velocity(Ntime,:)/10^3, 'LineWidth', P.LineWidth); grid on;
    if P.pxline ==1
    if i.grid.i_xpoint(1) > 1;  xline(X(i.grid.i_xpoint(1)),'-'); end
    if i.grid.i_xpoint(2) > 1;  xline(X(i.grid.i_xpoint(2)),'-'); end
    end
    if isstruct(P.solps);   hold on;
        %           species dependent(D0, D+1, C0, C+1:6)
        if P.avg == 0
            plot(solps.geom.dspar,-1*solps.datat.ua(:,2)/10^3,'--k', 'LineWidth', P.LineWidth)
        else
            plot(solps.geom.dspar(:,1),solps.data.ua(:,1,2)/10^3,'--k', 'LineWidth', P.LineWidth)
            tmpX = [solps.geom.dspar(1,1); solps.geom.dspar(:,1); flip(solps.geom.dspar(:,1))];
            tmpY = [solps.data.ua(1,1,2); solps.data.ua(:,2,2);  flip(solps.data.ua(:,3,2))]/10^3;
            patch(tmpX,tmpY,[0 0 0],'FaceAlpha',0.3,'LineStyle','none') % the errorbar
        end
    end
    if  keeplabel; vlab = ylabel('$v_{\mathrm{\|}}$ (10$^3$ms$^{-1}$)','FontSize',P.FontSize,'interpreter','latex');  end
    if max(P.ylim_v) ~=0;    ylim(P.ylim_v); end
    if max(P.ytic_v) ~=0;  yticks(P.ytic_v); end
    
if P.flip ==1
    if P.figtight ==0;    xlabel('distance to target (m)','FontSize',P.FontSize,'interpreter','latex'); % ,'interpreter','latex'
    else
        xticklabels({''}); xticks(P.xticks);
    end
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream (m)','interpreter','latex')
end
xlim(P.xlimits);

% neutral density
if P.figtight ==0
    figure(P.fignum+8);
    position = generalposition + [0 300 0 0 ];
    set(gcf, 'OuterPosition',position);
else
   if ~isfig; ax.nd = axes('Position',ndpos, 'FontSize', P.FontSize-1); 
   else; ax.nd = fig.Children(3+plimp);  ax.nd;  fig.CurrentAxes = ax.nd; end; box on;
end
if P.hold == 1;  hold on; end
plot(X,o.neutral_density(Ntime,:), 'LineWidth', P.LineWidth); grid on; 

if isstruct(P.solps);   hold on;
    if P.avg ==0
        plot(solps.geom.dspar(:,1),solps.data.na_correct,'--k', 'LineWidth', P.LineWidth)
    else
        plot(solps.geom.dspar(:,1),solps.data.na_correct(:,1),'--k', 'LineWidth', P.LineWidth)
        tmpX = [solps.geom.dspar(1,1); solps.geom.dspar(:,1); flip(solps.geom.dspar(:,1))];
        tmpY = [solps.data.na_correct(1,1); solps.data.na_correct(:,2);  flip(solps.data.na_correct(:,3))];
        patch(tmpX,tmpY,[0 0 0],'FaceAlpha',0.3,'LineStyle','none') % the errorbar
    end
end
if  keeplabel; ndlab = ylabel('$n_{\mathrm{D0}}$ (m$^{-3}$)','FontSize',P.FontSize,'interpreter','latex'); end
set(gca,'YScale','log')
if max(i.physics.initial_nb) > 0 % plot neutral background
    if P.plotbackground ==1;  hold on;
    if i.grid.i_xpoint(1) > 1
        plot(X([1,i.grid.i_xpoint(1)]),[1 1]*i.physics.initial_nb(1),...
            '--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth);
        plot(X([1,i.grid.i_xpoint(1)]),[1 1]*i.physics.initial_nb(2),...
            '--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth);
    end
    if i.grid.i_xpoint(2) > 1
        plot(X([i.grid.i_xpoint(1),i.grid.i_xpoint(2)]),...
            [1 1]*i.physics.initial_nb(3),...
            '--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth);
        if P.pxline == 1
        xline(X(i.grid.i_xpoint(2)), '-',  'X-point','FontSize',P.FontSize-2, 'interpreter','latex',...
            'LabelHorizontalAlignment', 'left')
        end
    end
    
     plot(X([i.grid.i_xpoint(2),end]),...
            [1 1]*i.physics.initial_nb(4),'--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth)
     leg2(1) = plot(X([i.grid.i_xpoint(2),end]),...
      [1 1]*i.physics.initial_nb(5),...
            '--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth);
      %legend(leg2(1),'$n_{\mathrm{b}}$','FontSize',P.FontSize,'interpreter','latex','location','best'); 
      hold off;
    end
end
if P.flip ==1
    if P.figtight ==0;    xlabel('distance to target (m)','FontSize',P.FontSize,'interpreter','latex'); % ,'interpreter','latex'
    else
        xticklabels({''}); xticks(P.xticks);
    end
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream (m)','interpreter','latex')
end
xlim(P.xlimits);
if P.ylim_a ~=0;    ylim(P.ylim_a); end
if P.ytic_a ~=0;  yticks(P.ytic_a); end


% neutral velocity
    if P.figtight ==0
        figure(P.fignum+7);
        position = generalposition + [0 600 0 0 ];
        set(gcf, 'OuterPosition',position);
    else
       if ~isfig; ax.nv = axes('Position',nvpos, 'FontSize', P.FontSize-1); else; ax.nv = fig.Children(2+plimp);  fig.CurrentAxes = ax.nv; end; ax.nv; box on;
    end
    if P.hold == 1;  hold on; end
    plot(X,o.neutral_velocity(Ntime,:)/10^3, 'LineWidth', P.LineWidth); grid on;
    if P.pxline ==1
    if i.grid.i_xpoint(1) > 1;  xline(X(i.grid.i_xpoint(1)),'-'); end
    if i.grid.i_xpoint(2) > 1;  xline(X(i.grid.i_xpoint(2)),'-'); end
    end
    if isstruct(P.solps);   hold on;
        %           species dependent(D0, D+1, C0, C+1:6)
        if P.avg == 0
            plot(solps.geom.dspar,-1*solps.datat.ua(:,1)/10^3,'--k', 'LineWidth', P.LineWidth)
        else
            plot(solps.geom.dspar(:,1),solps.data.ua(:,1,1)/10^3,'--k', 'LineWidth', P.LineWidth)
            tmpX = [solps.geom.dspar(1,1); solps.geom.dspar(:,1); flip(solps.geom.dspar(:,1))];
            tmpY = [solps.data.ua(1,1,1); solps.data.ua(:,2,1);  flip(solps.data.ua(:,3,1))]/10^3;
            patch(tmpX,tmpY,[0 0 0],'FaceAlpha',0.3,'LineStyle','none') % the errorbar
        end
    end
    if  keeplabel; nvlab = ylabel('$v_{\mathrm{\|D}}$ (10$^3$ms$^{-1}$)','FontSize',P.FontSize,'interpreter','latex');  end
    if max(P.ylim_nv) ~=0;    ylim(P.ylim_nv); end
    if max(P.ytic_nv) ~=0;  yticks(P.ytic_nv); end
    
if P.flip ==1
    if P.figtight ==0;    xlabel('distance to target (m)','FontSize',P.FontSize,'interpreter','latex'); % ,'interpreter','latex'
    else
        xticklabels({''}); xticks(P.xticks);
    end
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream (m)','interpreter','latex')
end
xlim(P.xlimits);

% molecules
    if P.figtight ==0
        figure(P.fignum+7);
        position = generalposition + [0 600 0 0 ];
        set(gcf, 'OuterPosition',position);
    else
       if ~isfig; ax.m = axes('Position',mpos, 'FontSize', P.FontSize-1); else; ax.m = fig.Children(1+plimp);  fig.CurrentAxes = ax.m; end; ax.m; box on;
    end
    if P.hold == 1;  hold on; end
    plot(X,o.molecule(Ntime,:), 'LineWidth', P.LineWidth); grid on;
    if P.pxline ==1
    if i.grid.i_xpoint(1) > 1;  xline(X(i.grid.i_xpoint(1)),'-'); end
    if i.grid.i_xpoint(2) > 1;  xline(X(i.grid.i_xpoint(2)),'-'); end
    end
    if isstruct(P.solps);   hold on;
        %           species dependent(D0, D+1, C0, C+1:6)
        if P.avg == 0
            plot(solps.geom.dspar,-1*solps.datat.dmb2(:,1),'--k', 'LineWidth', P.LineWidth)
        else
            plot(solps.geom.dspar(:,1),solps.data.dmb2(:,1,1),'--k', 'LineWidth', P.LineWidth)
            tmpX = [solps.geom.dspar(1,1); solps.geom.dspar(:,1); flip(solps.geom.dspar(:,1))];
            tmpY = [solps.data.dmb2(1,1,1); solps.data.dmb2(:,2,1);  flip(solps.data.dmb2(:,3,1))];
            patch(tmpX,tmpY,[0 0 0],'FaceAlpha',0.3,'LineStyle','none') % the errorbar
        end
    end
if max(i.physics.initial_mb) > 0; hold on; % plot neutral background;
    if P.plotbackground ==1 
    if i.grid.i_xpoint(1) > 1
        plot(X([1,i.grid.i_xpoint(1)]),[1 1]*i.physics.initial_mb(1),...
            '--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth);
        plot(X([1,i.grid.i_xpoint(1)]),[1 1]*i.physics.initial_mb(2),...
            '--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth);
    end
    if i.grid.i_xpoint(2) > 1
        plot(X([i.grid.i_xpoint(1),i.grid.i_xpoint(2)]),...
            [1 1]*i.physics.initial_mb(3),...
            '--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth);
        if P.pxline == 1
        xline(X(i.grid.i_xpoint(2)), '-',  'X-point','FontSize',P.FontSize-2, 'interpreter','latex',...
            'LabelHorizontalAlignment', 'left')
        end
    end   
     plot(X([i.grid.i_xpoint(2),end]),...
            [1 1]*i.physics.initial_mb(4),'--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth)
     leg(1) = plot(X([i.grid.i_xpoint(2),end]),...
      [1 1]*i.physics.initial_mb(5),...
            '--','color',[0, 0.4470, 0.7410],'LineWidth',P.LineWidth);
      %legend(leg(1),'$n_{\mathrm{b,D2}}$','FontSize',P.FontSize,'interpreter','latex','location','best'); 
    end
    hold off;
end
    if  keeplabel; mlab = ylabel('$n_{\mathrm{D2}}$ ($m^{-3}$)','FontSize',P.FontSize,'interpreter','latex');  end
    set(gca,'YScale','log')
    if max(P.ylim_m) ~=0;    ylim(P.ylim_m); end
    if max(P.ytic_m) ~=0;  yticks(P.ytic_m); end

if P.flip ==1
     xlabel('distance to target (m)','FontSize',P.FontSize,'interpreter','latex'); % ,'interpreter','latex'
     xticks(P.xticks);    
    set(gca,'XDir','reverse');
else
    xlabel('distance from upstream (m)','interpreter','latex')
end
xlim(P.xlimits);

if P.figtight ==1
    %            xp = 5.6; %5.82; %6;
    
    if  keeplabel
    ndlab.Position(1) =P.xp;
    qlab.Position(1) =P.xp;
    tlab.Position(1) = P.xp;
    nlab.Position(1) = P.xp;
    vlab.Position(1) =P.xp;
    nvlab.Position(1) =P.xp;
    mlab.Position(1) =P.xp;
    end
    pause(0.5)
end

end


% end
% P.save =1;
% if P.save ==1
%     savenames = {'q','te','ne','ui','na','fc'};
%     inccar = 5;
%     if isstruct(P.solps) inccar = 6; end
%     for i = 1:inccar
%         try
%             pause(1.0)
%     h = tightfig(figure(P.fignum+3+i));
%     tmpstr = strcat('Figures/',P.title,'_DvsS_R_',savenames{i});
%     saveas(h,tmpstr,'pdf');
%         catch
%     disp(' error 873')
%         end
%     end
% end
