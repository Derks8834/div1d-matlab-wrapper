% script to run a series of div1d runs from matlab

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023


clearvars; close all; clc;
addpath('div1d-matlab-wrapper');
addpath('div1d');
addpath('div1d/obj');


div1d_restart_file_path = '/home/emc/derks/Desktop/div1d/runs/TCV150678_v1.0.0_stat_avg_rad_vis5_neu1.8_qu23.6_gv4_R0.99_fR0.5_tn1_RL0.5/';

%% Compare with SOLPS-ITER: 
% This requires the div1d-matlab-toolbox 
% SOLPS shot(s)
shot_num = 150678;

% Method of SOLPS data extraction for DIV1D input.txt
%   maxqt   = flux tube through max q at target
%   maxqm   = flux tube trhough max q at midplane
%   avg     = average of tubes between above
%   maxq    = maximum of q everywhere along leg (legacy)
%   avg_rad = averaged with radial losses
solps_extr_meth = {'avg_rad'};% maxqt

%% General
user = 'derks';
date_str = date;  
tag = 'v2.0.0b';  % version of DIV1D (use git tag command to find it)
% versionX.Y.Z 
% X when backwards compatibility is lost. (when the namelist changes)
% Y when new features are added with compatibility.
% Z when a new version only contains small changes. 
% alpha addition when this version was under development
%% Settings of simulation/experiment
time = 10e-3;   % [s]
delta_t = 1e-6; % [s]
ntime = round(time/delta_t,0);
out_delta_t = 1e-4; % [s]
nout = round(out_delta_t/delta_t,0); % write output every nout steps
type = 'dyn'; % static runs to fit on with restart to see viscocity

puffing_source = 4*10^20;
% Recycling (generally higher for D than N)
base_R = 0.99;
% Fraction Redistributed
base_fR = 0.5;
% Tau 
base_tn = 1;
% arteficial viscosity
% viscosity = [0.04 0.2 0.6 1 3 5];
viscosity = 5;
% radial loss factor
RL = 0.5;

%% Loop
for imet = 1:length(solps_extr_meth)


method      = solps_extr_meth{imet};
solps_num   = shot_num;
shot_data_solps = strcat('./../SOLPS/data/TCV',num2str(solps_num));

% load solps data
load(strcat(shot_data_solps,'.mat')); 
tmp = SOLPS2DIV1D_dev(shot);


for iscan = 1 %:length(viscosity)
    R  = base_R;%list_R(iscan);% 0.8;
    fR = base_fR;%list_fR(iscan);%base_fR; %0.8; 
    tn = base_tn;%list_tn(iscan);%base_tn;%1d-4;
    eta = 5;%viscosity(iscan);
    
    try
        car_frc = list_car_frc(iscan);
    catch
        car_frc = round(median(tmp.(method).data.cratio(:,1)),3); 
    end
    %       data=read_namelist('input.test','INDATA');
%       fid=fopen('input.test_out','w');
%       write_fortran_namelist(fid,data,'INDATA');
%       fclose(fid);

% Load default input parameters
     [numpar,phypar] = default_input;

% Change  parameters
    numpar{01,1} = 'restart';                   numpar{01,2} = '.true.';% was: .true.
    numpar{02,1} = 'Nx';                        numpar{02,2} = 500; %default 500
    numpar{04,1} = 'ntime';                     numpar{04,2} = ntime;
    numpar{05,1} = 'nout';                      numpar{05,2} = nout;
    numpar{06,1} = 'delta_t';                   numpar{06,2} = delta_t;
    numpar{12,1} = 'viscosity';                 numpar{12,2} = eta; %default in code = 0.04
    
    phypar{01,1} = 'initial_n';                 phypar{01,2} = abs(tmp.(method).data.ne(1,1));%2.50*1e19;
    phypar{02,1} = 'initial_v';                 phypar{02,2} = 0.0d+5;
    phypar{03,1} = 'initial_T';                 phypar{03,2} = mean(tmp.(method).data.te(:,1));%10.4d+01
    phypar{04,1} = 'initial_a';                 phypar{04,2} = mean(tmp.(method).data.na_correct(:,1));%1d+13;
    phypar{05,1} = 'gamma';                     phypar{05,2} = 6.0d+0;
    phypar{06,1} = 'mass';                      phypar{06,2} = 3.3436d-27;
    phypar{07,1} = 'q_parX';                    phypar{07,2} = abs(tmp.(method).data.qpar(1,1));%2.935d+07;
    phypar{08,1} = 'Gamma_X';                   phypar{08,2} = 0.0d+24;
    phypar{09,1} = 'energy_loss_ion';           phypar{09,2} = 3.0d+1;
    phypar{10,1} = 'neutral_residence_time';    phypar{10,2} = tn;%1.0d-4;%         default for my runs: 1.0d-41d+20; % was 1.0d-4 in slide acht van pptx 1-28-2020 (DIV1D Density scan)
    phypar{11,1} = 'redistributed_fraction';    phypar{11,2} = fR;%0.8; %           default for my runs: 0.8 was 0.8 in the above
    phypar{12,1} = 'carbon_concentration';      phypar{12,2} = car_frc;% median less affected by high target 1d-2;
    phypar{13,1} = 'L';                         phypar{13,2} = round(tmp.(method).geom.dspar(1,1),2); %7.065;
    phypar{14,1} = 'recycling';                 phypar{14,2} = R;%74/100;%1; % default for my runs: 74/100;
    phypar{15,1} = 'sintheta';                  phypar{15,2} = round(tmp.(method).data.sintheta(1),3);%6.460d-02;
    phypar{16,1} = 'case_AMJUEL';               phypar{16,2} = '.true.';
    phypar{17,1} = 'minimum_temperature';       phypar{17,2} = 0.1d-0;
    phypar{18,1} = 'gas_puff_source';           phypar{18,2} = puffing_source;
    phypar{19,1} = 'gas_puff_location';         phypar{19,2} = 5.5d+0;
    phypar{20,1} = 'gas_puff_width';            phypar{20,2} = 1.0d+20;
    phypar{21,1} = 'elm_start_time';            phypar{21,2} = 1000;
    phypar{22,1} = 'elm_ramp_time';             phypar{22,2} = 500;
    phypar{23,1} = 'elm_time_between';          phypar{23,2} = 2d+3;
    phypar{24,1} = 'elm_expelled_heat';         phypar{24,2} = 0d+3;
    phypar{25,1} = 'elm_expelled_particles';    phypar{25,2} = 0d+18;
    phypar{26,1} = 'switch_elm_heat_flux';      phypar{26,2} = 0;
    phypar{27,1} = 'switch_elm_density';        phypar{27,2} = 0;
    phypar{28,1} = 'switch_elm_series';         phypar{28,2} = 0;
    phypar{29,1} = 'radial_loss_factor';        phypar{29,2} = RL;
%   phypar{30,1} = 'radial_loss_gaussian';      phypar{30,2} = 0;
%   phypar{31,1} = 'radial_loss_width';         phypar{31,2} = 7/6;
%   phypar{32,1} = 'radial_loss_location';      phypar{32,2} = 3.0;
    phypar{33,1} = 'switch_dyn_nu';             phypar{33,2} = 1;
    phypar{34,1} = 'switch_dyn_gas';            phypar{34,2} = 1;
    phypar{35,1} = 'switch_dyn_rec';            phypar{35,2} = 1;
    phypar{36,1} = 'switch_dyn_rad_los';        phypar{36,2} = 1;
    phypar{37,1} = 'switch_car_con_prf';        phypar{37,2} = 1;
    
  
    dyn_nu = ones(1,ntime)*phypar{01,2} + [zeros(1,floor(ntime/5)) ones(1,ntime-floor(ntime/5))]*1*10^18; 
    dyn_gas = ones(1,ntime)*phypar{18,2} + [zeros(1,floor(ntime/4)) ones(1,ntime-floor(ntime/4))]*1*10^20;
    dyn_rec = ones(1,ntime)*phypar{14,2} - [zeros(1,floor(ntime/2)) ones(1,ntime-floor(ntime/2))]*0.05 ;
    dyn_rad_los =ones(1,ntime)*phypar{29,2} - [zeros(1,floor(ntime/3)) ones(1,ntime-floor(ntime/3))]*0.05;
    
% Write input
    description = 'test dynamic inputs and carbon profile';
    write_input(numpar,phypar,'nu',dyn_nu,'gas',dyn_gas,'R',dyn_rec,'RL',dyn_rad_los,'description',description);
   
% Save logic
% directory logic
%     shot   vers.  type settings
% 'TCV150683_v1.0.0_stat_maxqt_vis0.9_nu2_qu30_gv0_R74_fRd80/'
%
% struct logic (.mat)
%  shot.type.settings. (information)
%  (information: tag, user, date, output)
    
% Define directory
% (e.g.'TCV150683_v1.0.0_stat_maxqt_vis0.9_nu2_qu30_gv0_R74_fRd80/')
%     strucdir = [method,'_',...
%             'vis', num2str(numpar{12,2}),'_', ... viscocity
%             'neu', num2str(round(phypar{01,2}/10^19,1)),'_',... ne upstream
%             'qu', num2str(round(phypar{07,2}/10^6,1)),'_',... qpar upstream
%             'gv', num2str(phypar{18,2}),'_',... Gas Valve
%             'R', num2str(phypar{14,2}),'_',...  Recycling
%             'fR', num2str(phypar{11,2}),'_',... fraction redistributed
%             'tn', num2str(phypar{10,2})];%  neutral residence time   
        
%     dir = ['TCV',num2str(solps_num),'_',... 
%             tag,'_',...
%             type,'_',...
%             strucdir];
%         dir = 'test_dyn_input_29_06';
% %     strucdirm = strrep(strucdir, '.', ''); % matlab structure field name
%         strucdirm = 'test_dyn_input_29_06';
        
% Run DIV1D
%     runnit = 0; % you can either run this directory or read the output, not yet sequentially. 
%     if runnit ==1
%         %run_program(dir,'run_div1d',1,'nohup',1); % it is possible to run on the background
%         run_program(dir,'run_div1d',1,'restart_file_path', div1d_restart_file_path );% just make the input files
%     else % read the output and save to .mat struct
%         output =  read_output(dir,'branch_path',div1d_branch_path);
%         savepath    = strcat('./DIV1D/',shot_data_solps,'.mat');
%         try 
%         load(savepath);
%         catch
%         disp('making new .mat')
%         end
%          shot.div1d.(type).(strucdirm).output   = output;
%          shot.div1d.(type).(strucdirm).tag      = tag;
%          shot.div1d.(type).(strucdirm).date     = date_str;
%          shot.div1d.(type).(strucdirm).usr      = user;
% 
%         save(savepath,'shot')
%     end
end    
end
 
