clear all;
clc;
close all; 
% add the div1d base folder to path
addpath(genpath('./div1d_astra'))
addpath(genpath('./'))

[output,input]= div1dread_v600('/home/unix/derks/Desktop/codes/div1d_public/div1d_astra/int-tests/12_TCV_dynamics/div1d_output.txt');


plotdiv1d_v600(output,input)

figure(88)
plot_div1d_polgeom(input)