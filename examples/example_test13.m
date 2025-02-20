clear all;
clc;
close all; 
% add the div1d base folder to path
addpath(genpath('./div1d_astra'))
addpath(genpath('./'))

path_div1d = '/home/unix/derks/Desktop/projects/dynamics/models/div1d-wall/div1d';

path12 = [path_div1d,'/int-tests/12_TCV_dynamics/div1d_output.txt'];
path13 = [path_div1d,'/int-tests/13_TCV_dyn_wall_association/div1d_output.txt'];

[output,input]= div1dread_v600(path12);
[output13,input13] = div1dread_v600(path13);

plotdiv1d_v600(output,input)

plotdiv1d_v600(output13,input13)

figure(88)
plot_div1d_polgeom(input)