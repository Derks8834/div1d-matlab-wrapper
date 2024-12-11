% script to run div1d standalone from matlab

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% Dec 2024


%% read one of the tests:

% first you have to compile the code and run the test
% make sure you have ifortran version 19 or later installed
% go into the div1d repository and use following commands
% Make ( if this produces an error, let me know)
%./test_div1d.sh
% Type 12: -> test 12 should run and take less than 5 minutes. 
% If now in int-test/12_TCV_dynamics/ you see div1d_output.txt, we are in business. 

% add div1d to path as the test is located there.
addpath(genpath('./div1d_astra'))
% current-directory; div1d-matlab-wrapper
addpath(genpath('./')) 

% read and plot the outcome
[output,input]= div1dread_v600('/home/unix/derks/Desktop/codes/div1d_public/div1d_astra/int-tests/12_TCV_dynamics/div1d_output.txt');
plotdiv1d_v600(output,input)
figure(88)
plot_div1d_polgeom(input)%% Settings of simulation/experiment

% post-process the outcome and perform code-checks
[proc] = process_div1d_output_v600(output,input,'ploterror',1,'plotbal',1);
% note that new versions of matlab support function arguments by default,
% however, in all repos here, we have our own varargin function argument
% parsing that also works on computer clusters which are not on matlab 2022
% or later.



%% setup a simulation of your own.

%%% Go out of the div1d-matlab-wrapper and have as sub-folders:
% div1d_astra
% div1d-matlab-wrapper
% then continue

run aug_simple_ppcf2024  % this is a simple case similar to the article

% setup what type of run you want
% do not run with core and neutral reservoirs
numpar.evolve_core = 0;
numpar.evolve_background = [0 0 0 0 0];
% do not run with neutral momentum and molecules 
numpar.evolve_neutral_momentum = 0;
numpar.evolve_molecule = 0.0;
% (div1d did not have these features when ppcf2024 was published)

% the following should be set to run a global particle balance simulation
% phypar.extern_molecule_density = ..
% phypar.core_confinement_time = 6.736841689E-02;
% phypar.core_ext_neutral_pump = [0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 1.441199180E+03];
% phypar.core_fuelling = 1.000000000E+01;
% phypar.extern_neutral_ex = [0.000000000E+00 0.000000000E+00 7.379585459E+03];
% phypar.extern_molecule_ex = [0.000000000E+00 0.000000000E+00 5.797949912E+02];
% phypar.puff_rate_neutral =[0.000000000E+00 0.000000000E+00 -2.729660830E+20 -1.816175340E+21 3.489192769E+21];
% phypar.puff_rate_molecule =[0.000000000E+00 0.000000000E+00 6.207295765E+20 -1.483974546E+21 -1.751297196E+20];
% phypar.pump_rate_n =[0.000000000E+00 0.000000000E+00 1.000000000E+00 1.000000000E+00 1.000000000E+00];
% phypar.pump_rate_m =[0.000000000E+00 0.000000000E+00 1.000000000E+00 1.000000000E+00 1.000000000E+00];

% write your input files (these are going in your current directory)
[proceed] =write_input(numpar,phypar); 
% for dynamic simulations you have to specify time-varying inputs in dat
% files
% dat = struct; 
%       [proceed] =write_input(numpar,phypar,'dat',dat);
% type help write_input for more information (or look in the function)
% Without errors, proceed = true

% define simulation folder and directories
div1d_install_path = '/home/unix/derks/Desktop/codes/div1d_public/div1d_astra/'; % point to your instalation
simulation_path = './aug_simple_ppcf2024/';
runnit = 1; % turn to 0 to first check if all files are written correctly
nohup = 0; % to have div1d write into nohup.out on simulation progress
datinput = 0; % set to 1 to look for time-dependent .dat files and parse them to div1d
if proceed
run_div1d_code('','runs_path',simulation_path,"div1d_install_path",div1d_install_path,'run_div1d',runnit,'nohup',nohup,'datinput',datinput); 
end
% note that new versions of matlab support function arguments by default,
% however, in all repos here, we have our own varargin function argument
% parsing that also works on computer clusters which are not on matlab 2022
% or later.

%% wait a while
% input('wait for the simulation to finish'); 
% in this case nohup =0 so you will see when it is finished.

%% inspect outcome

[output,input]= div1dread_v600(strcat(simulation_path,'div1d_output.txt'));
plotdiv1d_v600(output,input)
%figure(88) % no geometry was properly set for this demo
%plot_div1d_polgeom(input)%% Settings of simulation/experiment

% post-process the outcome and perform code-checks
[proc] = process_div1d_output_v600(output,input,'ploterror',1,'plotbal',1);

%% now if you want to simulate with the reservoirs (!!WIP!!)
% obviously you want to keep in steady state the scrape-off layer solution
% in that case try:
[ ~,reservoir_settings] = get_div1d_chamber_params(input,output,...
              'core_fuelling',1e1,'core_density',phypar.initial_n*4,...
                                    'molecule_puff',[0 0 0 0.5 0.5]*1e22,...
                                    'atom_puff',[0 0 0 0 0],...
                                    'ato_pump',[0 0 1 1 1]*40,...
                                    'mol_pump',[0 0 1 1 1]*40,...
                                    'print',1); 

% this will print for you the settings that you need to obtain a stationary
% SOL in the vicinity of your original solution without evolving the
% reservoirs WARNING THIS MAY BE VERY UNPHYSICAL
% pasted the outcome
phypar.initial_ncore = 4.000000000E+19;
phypar.core_confinement_time = 2.222222222E-02;
phypar.core_ext_neutral_pump = [0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 2.078860825E+02];
phypar.core_fuelling = 1.000000000E+01;
phypar.extern_neutral_ex = [0.000000000E+00 0.000000000E+00 7.379585459E+03];
phypar.extern_molecule_ex = [0.000000000E+00 0.000000000E+00 5.797949912E+02];
phypar.puff_rate_neutral =[0.000000000E+00 0.000000000E+00 3.804314082E+21 0.000000000E+00 4.067216865E+23];
phypar.puff_rate_molecule =[0.000000000E+00 0.000000000E+00 9.856514850E+18 5.000000000E+21 5.000000000E+21];
phypar.pump_rate_n =[0.000000000E+00 0.000000000E+00 0.000000000E+00 4.963966973E+04 0.000000000E+00];
phypar.pump_rate_m =[0.000000000E+00 0.000000000E+00 0.000000000E+00 6.354460318E+03 1.803751804E+02];


% simulate again
numpar.evolve_core = 1; % this is a very bad core model...
numpar.evolve_background = [0 0 1 1 1];

simulation_path2 = './aug_simple_ppcf2024_reservoir/';

[proceed] =write_input(numpar,phypar); 
if proceed
run_div1d_code('','runs_path',simulation_path2,"div1d_install_path",div1d_install_path,'run_div1d',runnit,'nohup',nohup,'datinput',datinput); 
end
% this may take a while .... 


%%
[output2,input2]= div1dread_v600(strcat(simulation_path2,'div1d_output.txt'));

%% inspect the outcome
close all;
plotdiv1d_v600(output,input)
plotdiv1d_v600(output2,input2)
% figure(88) % no geometry was properly set for this demo
% plot_div1d_polgeom(input2)%% Settings of simulation/experiment

% post-process the outcome and perform code-checks
[proc2] = process_div1d_output_v600(output2,input2,'ploterror',1,'plotbal',1);
