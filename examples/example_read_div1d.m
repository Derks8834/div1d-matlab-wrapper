% script to read a series of div1d runs from matlab
% the standard case reads the unit-test runs
% go into the div1d and run ./test_div1d.sh to generate the outputs.

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

clearvars;  clc; 

% specify DIV1D install path
div1d_install_path   = '/home/unix/derks/Desktop/codes/div1d_neutrals/div1d/';
basepath =  [div1d_install_path,'unit-tests/'];

% add paths
addpath('functions_read','functions_plot','functions_analysis',basepath);

% select a div1d_output.txt from one of the tests
w_dir = pwd; 
cd(basepath)
[file,pathis] = uigetfile('*.txt');
cd(w_dir);
if file ~=0
    path = strcat(pathis,file);
end

% read the file
[output,input] = read_output(path);


% plot the output
plot_div1d(output,input);

%%
proc =  process_div1d_output(output,input);

%%
plot_div1d_error(output.X,proc.Derror)
