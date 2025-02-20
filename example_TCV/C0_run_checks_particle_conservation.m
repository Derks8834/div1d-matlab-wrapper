% script to run some tests 
% g.l.derks

% Fitting DIV1D from stagnation to target point for TCV with molecules
clearvars;
close all;
clc;

set(0, 'defaultFigureRenderer', 'painters')
 
%% get default input for DIV1D simulation
dir = struct;
device = 'TCV';
for check_iteration = 1:16
     div1d_install_path   = '/home/unix/derks/Desktop/projects/dynamics/models/div1d/';
%         div1d_install_path   = '< point to the folder where you cloned div1d > ';
        input_file_path      = pwd;
        restart_file_path    = '';
        branch_path          = '';
         runs_path            = '/home/unix/derks/Desktop/projects/dynamics/models/div1d/runs/gijsderks/NF2025/div1d/';
%         runs_path            = '< point the the folder where you want the simulations to be stored> ';
        dir.base = strcat(device);
        tmpstr = strcat('/',device,'_debug_check',num2str(check_iteration));
        dir.folders{check_iteration} = tmpstr;
end

%% run DIV1D
% use runnit = 1 to run but see if the folders are created correctly with runnit = 0
runnit = 1; 
for check_iteration =1:4 
            % get input
            run C0_inputs_check_particle_conservation_sol  
      
            % writes input.txt and possibly .dat files      
            nohup = 1;         
            format long;
            write_input(numpar,phypar) %,'dat',dat);

            % run the code
            run_div1d_code(strcat(dir.base,dir.folders{check_iteration}),"runs_path",runs_path,"div1d_install_path",div1d_install_path,'run_div1d',runnit,'nohup',nohup,'datinput',datinput);       
            
end
test = input('wait for the simulation to finish');
%% Look at output
for check_iteration = 1:4 
        [out{check_iteration},indiv{check_iteration}] =read_output('rnd','runs_path',[pwd,'/data/div1d/',dir.base, dir.folders{check_iteration}],'trymatfile',true);
        [proc{check_iteration}] = process_div1d_output_v600(out{check_iteration},indiv{check_iteration},'ploterror',0,'plotbal',1);
                 
        plotdiv1d = true;
        if plotdiv1d; plotdiv1d_v600(out{check_iteration},indiv{check_iteration},'hold', 1,'fignum',30); end

  test = input('look at next simulation');
end
