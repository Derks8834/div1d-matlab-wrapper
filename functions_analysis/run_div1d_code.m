function run_div1d_code(dir, varargin)%input_file_path,restart_file_path,branch_path,dir_name)
% run_div1d(dir_name,varargin)
% run_div1d(dir_name,'div1d_install_path','/home/emc/derks/Desktop/div1d/')
% SHORT:
%  runs DIV1D in a given directory 'dir'
% INPUT:        
%  dir      == string folder "dir" is made to run and store data in
%  varargin == 'div1d_install_path','string' (def == '/home/emc/derks/Desktop/div1d/'') 
%           == 'input_file_path','string'    (def == '')
%           == 'restart_file_path','string'  (def == '')
%           == 'runs_path','string'          (def == [div1d_install_path,'runs/']) (path where runs are stored)
%           == 'run_div1d'                   (def == 0 ) (NOTE:carefull with long runs)
%           == 'nohup'                       (def == 0) (NOTE: carefull!)
%           == 'datinput'                    (def == 0) (input ".dat" files)
%
%  NOTE: data is stored in directory: [runs_path,dir]
%
% OUTPUT:
%   div1d_output.txt
%   div1d_restart_new.txt
%   [files written in dir]
% use read_output.m to read the .txt files
% use process_div1d_output.m to interpret the raw output.

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

%% Handle parameters
D.div1d_install_path   = '/home/emc/derks/Desktop/div1d/';
D.input_file_path      = './'; % path where DIV1D is installed
D.restart_file_path    = ''; % path with restart file
D.runs_path             = ''; % path where "dir" folder is stored
D.run_div1d            = 0;
D.datinput             = 0; % vector inputs
D.nohup                = 0; % runs in background CAREFULL
P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end
if isempty(P.runs_path) % if no runs folder path is given explicitly
    P.runs_path = [P.div1d_install_path,'runs/']; % run from div1d dir
end
%% Construct commands

% make directory to run div1d
command_0 = ['bash -c "mkdir -p ',P.runs_path,dir,'"'];
% copy input.txt file 
command_1 = ['bash -c "cp ',P.input_file_path,'input.txt ',...
    P.runs_path,dir,'/input.txt','"'];
command_1A = ['bash -c "mv ',P.input_file_path,'/*.dat ',...
    P.runs_path,dir,'/','"'];
% copy div1d_restart_old.txt file
if ~isempty(P.restart_file_path)
    command_2 = ['bash -c "cp ', P.restart_file_path,'div1d_restart_new.txt ',...
               P.runs_path,dir,'/div1d_restart_old.txt','"'];    
end

% run div1d
if P.nohup == 1
    command_3 = ['bash -c "cd ',P.runs_path,dir,...
                 '; rm nohup.out; rm div1d_output.mat;' ,...
                 'nohup ',P.div1d_install_path,'obj/div1d.exe < input.txt & "'];
else
    command_3 = ['bash -c "cd ',P.runs_path,dir,...
                '; ',P.div1d_install_path,'obj/div1d.exe < input.txt ; exit"'];
end
%% execute commmands
system(command_0);
system(command_1);
if P.datinput ==1
    system(command_1A);
end
if ~isempty(P.restart_file_path)
    system(command_2);
end
if P.run_div1d ==1
    system(command_3);
end
end
