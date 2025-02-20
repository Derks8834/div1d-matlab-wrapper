% TCV script to run the div1d stationary comparison with SOLPS-ITER
% g.l.derks

% Fitting DIV1D from stagnation to target point for TCV with molecules
clearvars;
close all;
clc;
% addsearchpaths
set(0, 'defaultFigureRenderer', 'painters')

basepath = '/home/unix/derks/Desktop/projects/dynamics';
datpath = [basepath '/data/'];
e_charge = 1.61e-19;
%% TCV SOLPS simulations
plot_solps = false;
run A11_get_TCV_solps_information % gives .dat files

%% get default input for DIV1D simulation
dir = struct;

for check_iteration = 1:14
        div1d_install_path   = '/home/unix/derks/Desktop/projects/dynamics/models/div1d-debug/div1d/';
        input_file_path      = pwd;
        restart_file_path    = '';
        branch_path          = '';
        runs_path            = '/home/unix/derks/Desktop/projects/dynamics/models/div1d/runs/gijsderks/NF2025/div1d/';
        dir.base = strcat(device);
        tmpstr = strcat('/',device,'_debug_check',num2str(check_iteration));
        dir.folders{check_iteration} = tmpstr;
end

%% run DIV1D
% this you can set in get_tcv_default_div1d_input.m to push the steady state.
runnit = 0;
for check_iteration = 1:14 %13:14 %10:12  %1:8 %:5 %1:3 %1:6
            run C2_get_TCV_default_div1d_input_checks     
            % writes input.txt and possibly .dat files      
            nohup = 1;         
            format long;
            write_input(numpar,phypar,'dat',dat);
            run_div1d_code(strcat(dir.base,dir.folders{check_iteration}),"runs_path",runs_path,"div1d_install_path",div1d_install_path,'run_div1d',runnit,'nohup',nohup,'datinput',datinput);       
            test = input('wait for the simulation to finish');
end

%%
% read_dir = '' write directories manually for the different tests

%solps_select = 4;
for check_iteration =  8:14%10:12 %5:6
    %:6  % this is the one where we test settings
        % look at output and produce follow-up set
        [out{check_iteration},indiv{check_iteration}] =read_output('rnd','runs_path',[pwd,'/data/div1d/',dir.base, dir.folders{check_iteration}],'trymatfile',true);
        [proc{check_iteration}] = process_div1d_output_v600(out{check_iteration},indiv{check_iteration},'ploterror',0,'plotbal',0);
                 
        % given outcome in out{1} and indiv{1} -> calculate chamber params
        % with settings for core fuelling, density, molecule puff and atom puff
%         [~,set_] = get_div1d_chamber_params(indiv{check_iteration},out{check_iteration},...
%                                     'core_fuelling',1e21,'core_density',indiv{check_iteration}.physics.initial_ncore,...
%                                     'molecule_puff',[0 0 0 0.5 0.5]*inputs.puff(solps_select),...
%                                     'atom_puff',[0 0 0 0 0], 'print',1); 
        plotdiv1d = false;% true;
        if plotdiv1d; plotdiv1d_v600(out{check_iteration},indiv{check_iteration},'hold', 1,'fignum',30); end


  test = input('look at next simulation');
end
%%
% the one dynamic simulation (requires
check_iteration = 3;
[out{check_iteration},indiv{check_iteration}] =read_output('rnd','runs_path',[pwd,'/data/div1d/',dir.base, dir.folders{check_iteration}],'trymatfile',true);
[proc{check_iteration}] = process_div1d_output_v600(out{check_iteration},indiv{check_iteration},'ploterror',1,'plotbal',1);
plotdiv1d = false;
if plotdiv1d; plotdiv1d_v600(out{check_iteration},indiv{check_iteration},'hold', 1,'fignum',30); end


%% FRF analysis
div1d = struct;
div1d.ms = MSstruct;
div1d.ms.t1 = 0.25;
div1d.ms.P = 1; % only 1 in analysis
div1d.ms.fs = MSstruct.fs/10;

o = out{9};
in = indiv{9}.dynamic;
% input
input.gas = frf_cutsignal(in.dyn_molecule_puff(:,4),o.time,div1d.ms.t1, div1d.ms.f0, div1d.ms.fs,div1d.ms.P);
% sol 
[div1d.T6] =frf_cutsignal(o.temperature(:,end-6)    ,o.time,div1d.ms.t1, div1d.ms.f0, div1d.ms.fs,div1d.ms.P);
% core
[div1d.nc] =frf_cutsignal(o.core_density(:)         ,o.time,div1d.ms.t1, div1d.ms.f0, div1d.ms.fs,div1d.ms.P);
% neutrals
[div1d.na4] =frf_cutsignal(o.extern_neutral_density(:,4),o.time,div1d.ms.t1, div1d.ms.f0, div1d.ms.fs,div1d.ms.P);
[div1d.na3] =frf_cutsignal(o.extern_neutral_density(:,3),o.time,div1d.ms.t1, div1d.ms.f0, div1d.ms.fs,div1d.ms.P);
terms = {'T6','nc','na4'};

%% DFT from input to outputs
FRF_D.T6 = frf_dft_main(input.gas,  div1d.T6,div1d.ms.f0,div1d.ms.f_lines,div1d.ms.P,0);
FRF_D.nc = frf_dft_main(input.gas,  div1d.nc,div1d.ms.f0,div1d.ms.f_lines,div1d.ms.P,0);
FRF_D.na4 = frf_dft_main(input.gas, div1d.na4,div1d.ms.f0,div1d.ms.f_lines,div1d.ms.P,0);
FRF_D.na3 = frf_dft_main(input.gas, div1d.na3,div1d.ms.f0,div1d.ms.f_lines,div1d.ms.P,0);

%% plot response
varlist =  {'T6','nc','na4','na3'};
colorlist = jet(4);
figure(101)
for ii=1:length(varlist)
    %magnitude
    subplot(211)
    if numpar.evolve_core==0
    title('frequency response without core')
else
    title('frequency response including core')
end
    f_exc = div1d.ms.f0*div1d.ms.f_lines; %f_excd; %  FRF_D.(varlist{ii}).f_exc;
    norm = abs(FRF_D.(varlist{ii}).TF(1));
    loglog(f_exc,abs(FRF_D.(varlist{ii}).TF)/norm,'o-','LineWidth',1); hold on;
    ylabel('magnitude')
    grid on;
    %phase
    subplot(212)
    Angle=unwrap(angle(FRF_D.(varlist{ii}).TF));
    Angle_deg=rad2deg(Angle);
    ylabel('phase')
    cor = 0;
    % sign correction for inverse relations
    %if ii == 1 && numpar.evolve_core ==0; cor = -180; end
    if ii == 1 && numpar.evolve_core ==1; cor = -180; end
    %if ii == 2 && numpar.evolve_core ==1; cor = -180; end
    %if ii == 4 && numpar.evolve_core ==1; cor = -180; end
    l(ii) =semilogx(f_exc,Angle_deg+cor,'o-','LineWidth',1); hold on;
    grid on;
    ylim([-180 0])
end

legend(l,varlist{:})
%%
savefigpath = '/home/unix/derks/Desktop/projects/dynamics/analysis/TCV/figures';
if numpar.evolve_core==0
    savestr = strcat(savefigpath,'/frf_wo_core_first');
else
    savestr = strcat(savefigpath,'/frf_w_core_first');
end

saveas(gcf,strcat(savestr,'.fig'));
saveas(gcf,strcat(savestr,'.png'));

%% look at the basis of the solutions
[Nt, Nx] = size(o.temperature);
insel = 6:Nt;

XT = transpose(o.temperature(insel,:));
Xn = transpose(o.density(insel,:));
Xq = transpose(o.q_parallel(insel,:));

[UT,ST,VT] = svd(XT);
[Un,Sn,Vn] = svd(Xn);
[Uq,Sq,Vq] = svd(Xq);
size(UT);

%% plot outcome
figure(77);
% for ii = 1:4
r = 3;
% subplot(3,1,1
plot(o.X,UT(:,1:r)); 
title('basis functions')
ylabel('temperature (eV)'); grid on;

% subplot(3,1,2)
% plot(o.X,Un(:,1:r)); 
% ylabel("n")
% subplot(3,1,3)
% plot(o.Xcb,Uq(:,1:r)); 
% ylabel("q")
xlabel('distance from upstream (m)')


hold off;
%%
rUT = transpose(UT(:,1:r));
rUn = transpose(Un(:,1:r));
rUq = transpose(Uq(:,1:r));

coT = rUT*XT; 
con = rUn*Xn;
coq = rUq*Xq;

[divpod.T1] =frf_cutsignal(coT(1,:)   ,o.time(insel),div1d.ms.t1, div1d.ms.f0, div1d.ms.fs,div1d.ms.P);
[divpod.T2] =frf_cutsignal(coT(2,:)   ,o.time(insel),div1d.ms.t1, div1d.ms.f0, div1d.ms.fs,div1d.ms.P);
[divpod.T3] =frf_cutsignal(coT(3,:)   ,o.time(insel),div1d.ms.t1, div1d.ms.f0, div1d.ms.fs,div1d.ms.P);
frfpod = struct;
frfpod.T1 = frf_dft_main(input.gas,  divpod.T1,div1d.ms.f0,div1d.ms.f_lines,div1d.ms.P,1);

%% Dynamic Mode Decomposition
x = XT(:,500:end-1);
xp = XT(:,500+1:end);
u = in.dyn_molecule_puff(insel,4);
u = [1e-20*u(500:end-1)']; % normalize and cut
r = 3;
p = r+1;

[A_tilde,B_tilde,U_hat] = DMDcinator(x,xp,u,r,p);

dep_prof = U_hat*B_tilde;
boptions = bodeoptions; boptions.FreqUnits = {'Hz'};
tf_tilde = tf(ss(A_tilde,B_tilde,eye(r),[],1e-4)); figure;bode(tf_tilde,boptions);

A_full = U_hat*A_tilde*pinv(U_hat); imagesc(A_full);

xx = x;
uu = u;
x0 = xx(:,1);
a0 = pinv(U_hat)*x0;

a = nan(r,size(xx,2));
a(:,1) = a0;
for i = 1 : size(a,2)
    a(:,i+1) = A_tilde*a(:,i)+B_tilde*uu(:,i);
end
x_est = U_hat*a;

makevideo = true;
videopath = '/home/unix/derks/Desktop/projects/dynamics/analysis/TCV/figures';
msg = '';
if isempty(msg);[dirstat,msg,msgID] = mkdir(strcat(videopath,dir.folders{9})); end
videofile = strcat(videopath,dir.folders{9},'/div1dvideo.gif');

figure(7)
i = 1;
h1 = plot(o.X,xx(:,i),'k','LineWidth',2);
hold('on');
h2 = plot(o.X,x_est(:,i),'--r','LineWidth',2);
xlabel('distance from upstream [m]');
ylabel('temperature [eV]');
hold off;
ylim([0,50])

gif(videofile,'DelayTime',1/12,'LoopCount',15);
for i = 1 : 10 : size(xx,2)
    h1.YData = xx(:,i);
    h2.YData = x_est(:,i);   
    drawnow;
    gif;
end
web(videofile) % view video
%%
figure(8)
plot(o.time(insel(1:end-500)),xx(400,:),'k'); hold on;
plot(o.time(insel(1:end-499)),x_est(400,:),'--r'); hold off;
legend('DIV1D','POD')
xlim([0.2 0.4]);
xlabel('time (s)');
ylabel('temperature in volume 400 (eV)')