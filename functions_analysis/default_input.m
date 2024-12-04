function [numpar,phypar] = default_input(varargin)
% Generates default input for DIV1D
% variable arguments:
%    - input_div1d('radial_losses',1,'time_dependent',1) (Default option)
%   allows to generate inputs according to previous versions of DIV1D
%   that did not have radial losses or time-dependent .dat inputs

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023


D.version = 'v4.0.1';
D.radial_losses = 1;
D.time_dependent= 1;
P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end
switch P.version
    case 'v1.0.0'
        
%% Numerical parameters ($div1d_numerics)
numpar{01,1} = 'restart';                   numpar{01,2} = '.true.';
numpar{02,1} = 'Nx';                        numpar{02,2} = 500;
numpar{03,1} = 'dxmin';                     numpar{03,2} = 0.05d-0;
numpar{04,1} = 'ntime';                     numpar{04,2} = 15000;
numpar{05,1} = 'nout';                      numpar{05,2} = 100;
numpar{06,1} = 'delta_t';                   numpar{06,2} = 1.0d-6;
numpar{07,1} = 'abstol';                    numpar{07,2} = 1.0d-6;
numpar{08,1} = 'reltol';                    numpar{08,2} = 1.0d-6;
numpar{09,1} = 'method';                	numpar{09,2} = 227;
numpar{10,1} = 'max_step';                  numpar{10,2} = 100000;
numpar{11,1} = 'max_attempts';              numpar{11,2} = 1000;
numpar{12,1} = 'viscosity';                 numpar{12,2} = 4.0d-2;
numpar{13,1} = 'evolve_density';            numpar{13,2} = 1;
numpar{14,1} = 'evolve_momentum';           numpar{14,2} = 1;
numpar{15,1} = 'evolve_energy';             numpar{15,2} = 1;
numpar{16,1} = 'evolve_neutral';            numpar{16,2} = 1;
numpar{17,1} = 'switch_density_source';     numpar{17,2} = 1.0d+0;
numpar{18,1} = 'switch_momentum_source';    numpar{18,2} = 1.0d+0;
numpar{19,1} = 'switch_energy_source';      numpar{19,2} = 1.0d+0;
numpar{20,1} = 'switch_neutral_source';     numpar{20,2} = 1.0d+0;
numpar{21,1} = 'switch_charge_exchange';    numpar{21,2} = 1.0d+0;
numpar{22,1} = 'switch_ionization';         numpar{22,2} = 1.0d+0;
numpar{23,1} = 'switch_excitation';         numpar{23,2} = 1.0d+0;
numpar{24,1} = 'switch_convective_heat';    numpar{24,2} = 1.0d+0;
numpar{25,1} = 'switch_recombination';      numpar{25,2} = 1.0d+0;
numpar{26,1} = 'switch_impurity_radiation'; numpar{26,2} = 1.0d+0;
numpar{27,1} = 'simple_sol';                numpar{27,2} = '.false.';
numpar{28,1} = 'filter_sources';            numpar{28,2} = '.false.';
%% Physics parameters ($div1d_physics)
phypar{01,1} = 'initial_n';                 phypar{01,2} = 3.1d+19;
phypar{02,1} = 'initial_v';                 phypar{02,2} = 0.0d+5;
phypar{03,1} = 'initial_T';                 phypar{03,2} = 1.044d+2;
phypar{04,1} = 'initial_a';                 phypar{04,2} = 1.0d+13;
phypar{05,1} = 'gamma';                     phypar{05,2} = 6.0d+0;
phypar{06,1} = 'mass';                      phypar{06,2} = 3.3436d-27;
phypar{07,1} = 'q_parX';                    phypar{07,2} = 0.5d+8;
phypar{08,1} = 'flux_expansion';            phypar{08,2} = 1.0d+0; 
%  B_field = 1/ ( 1 + (flux_expansion-1)*x/L)
phypar{09,1} = 'Gamma_X';                   phypar{09,2} = 0.0d+24;
phypar{10,1} = 'energy_loss_ion';           phypar{10,2} = 3.0d+1;
phypar{11,1} = 'neutral_residence_time';    phypar{11,2} = 1d+20;
phypar{12,1} = 'redistributed_fraction';    phypar{12,2} = 0.0d+0;
phypar{13,1} = 'carbon_concentration';      phypar{13,2} = 1d-2;
phypar{14,1} = 'L';                         phypar{14,2} = 2.0d+1;
phypar{15,1} = 'recycling';                 phypar{15,2} = 1.00d+0;
phypar{16,1} = 'sintheta';                  phypar{16,2} = 1.0d-1;
phypar{17,1} = 'case_AMJUEL';               phypar{17,2} = '.true.';
phypar{18,1} = 'minimum_temperature';       phypar{18,2} = 0.1d-0;
phypar{19,1} = 'gas_puff_source';           phypar{19,2} = 0.0d+24;
phypar{20,1} = 'gas_puff_location';         phypar{20,2} = 5.5d+0;
phypar{21,1} = 'gas_puff_width';            phypar{21,2} = 1.0d+20;
phypar{22,1} = 'elm_start_time';            phypar{22,2} = 1000;
phypar{23,1} = 'elm_ramp_time';             phypar{23,2} = 200;
phypar{24,1} = 'elm_time_between';          phypar{24,2} = 2d+3;
phypar{25,1} = 'elm_expelled_heat';         phypar{25,2} = 40d+3;
phypar{26,1} = 'elm_expelled_particles';    phypar{26,2} = 5d+18;
phypar{27,1} = 'switch_elm_heat_flux';      phypar{27,2} = 0;
phypar{28,1} = 'switch_elm_density';        phypar{28,2} = 0;
phypar{29,1} = 'switch_elm_series';         phypar{29,2} = 0;
if P.radial_losses ==1
    phypar{30,1} = 'radial_loss_factor';        phypar{30,2} = 0.0d+00;
    phypar{31,1} = 'radial_loss_gaussian';      phypar{31,2} = 0;
    phypar{32,1} = 'radial_loss_width';         phypar{32,2} = 1.0d+20;
    phypar{33,1} = 'radial_loss_location';      phypar{33,2} = 0.0d+00;
end
if P.radial_losses ==1
    if P.time_dependent ==1
    phypar{34,1} = 'switch_dyn_nu';             phypar{34,2} = 0;
    phypar{35,1} = 'switch_dyn_gas';            phypar{35,2} = 0;
    phypar{36,1} = 'switch_dyn_rec';            phypar{36,2} = 0;
    phypar{37,1} = 'switch_dyn_rad_los';        phypar{37,2} = 0;
    phypar{38,1} = 'switch_car_con_prf';        phypar{38,2} = 0;
    phypar{39,1} = 'switch_dyn_qpar';           phypar{39,2} = 0;
    phypar{40,1} = 'switch_dyn_red_frc';        phypar{40,2} = 0;
    end
end
       
    case 'v3.0.2'
        
    %% Numerical parameters ($div1d_numerics)
numpar{01,1} = 'restart';                   numpar{01,2} = '.true.';
numpar{02,1} = 'Nx';                        numpar{02,2} = 500;
numpar{03,1} = 'dxmin';                     numpar{03,2} = 0.05d-0;
numpar{04,1} = 'ntime';                     numpar{04,2} = 15000;
numpar{05,1} = 'nout';                      numpar{05,2} = 100;
numpar{06,1} = 'delta_t';                   numpar{06,2} = 1.0d-6;
numpar{07,1} = 'abstol';                    numpar{07,2} = 1.0d-6;
numpar{08,1} = 'reltol';                    numpar{08,2} = 1.0d-6;
numpar{09,1} = 'method';                	numpar{09,2} = 227;
numpar{10,1} = 'max_step';                  numpar{10,2} = 100000;
numpar{11,1} = 'max_attempts';              numpar{11,2} = 1000;
numpar{12,1} = 'viscosity';                 numpar{12,2} = 4.0d-2;
numpar{13,1} = 'evolve_density';            numpar{13,2} = 1;
numpar{14,1} = 'evolve_momentum';           numpar{14,2} = 1;
numpar{15,1} = 'evolve_energy';             numpar{15,2} = 1;
numpar{16,1} = 'evolve_neutral';            numpar{16,2} = 1;
numpar{17,1} = 'switch_density_source';     numpar{17,2} = 1.0d+0;
numpar{18,1} = 'switch_momentum_source';    numpar{18,2} = 1.0d+0;
numpar{19,1} = 'switch_energy_source';      numpar{19,2} = 1.0d+0;
numpar{20,1} = 'switch_neutral_source';     numpar{20,2} = 1.0d+0;
numpar{21,1} = 'switch_charge_exchange';    numpar{21,2} = 1.0d+0;
numpar{22,1} = 'switch_ionization';         numpar{22,2} = 1.0d+0;
numpar{23,1} = 'switch_excitation';         numpar{23,2} = 1.0d+0;
numpar{24,1} = 'switch_convective_heat';    numpar{24,2} = 1.0d+0;
numpar{25,1} = 'switch_recombination';      numpar{25,2} = 1.0d+0;
numpar{26,1} = 'switch_impurity_radiation'; numpar{26,2} = 1.0d+0;
numpar{27,1} = 'simple_sol';                numpar{27,2} = '.false.';
numpar{28,1} = 'filter_sources';            numpar{28,2} = '.false.';
%% Physics parameters ($div1d_physics)
phypar{01,1} = 'initial_n';                 phypar{01,2} = 3.1d+19;
phypar{02,1} = 'initial_v';                 phypar{02,2} = 0.0d+5;
phypar{03,1} = 'initial_T';                 phypar{03,2} = 1.044d+2;
phypar{04,1} = 'initial_a';                 phypar{04,2} = 1.0d+13;
phypar{05,1} = 'gamma';                     phypar{05,2} = 6.0d+0;
phypar{06,1} = 'mass';                      phypar{06,2} = 3.3436d-27;
phypar{07,1} = 'q_parX';                    phypar{07,2} = 0.5d+8;
phypar{08,1} = 'flux_expansion';            phypar{08,2} = 1.0d+0; 
%  B_field = 1/ ( 1 + (flux_expansion-1)*x/L)
phypar{09,1} = 'Gamma_X';                   phypar{09,2} = 0.0d+24;
phypar{10,1} = 'energy_loss_ion';           phypar{10,2} = 3.0d+1;
phypar{11,1} = 'neutral_residence_time';    phypar{11,2} = 1d+20;
phypar{12,1} = 'redistributed_fraction';    phypar{12,2} = 0.0d+0;
phypar{13,1} = 'num_impurities';            phypar{13,2} = 5;
phypar{14,1} = 'impurity_concentration';    phypar{14,2} = [0.01 0 0 0 0];
phypar{15,1} = 'impurity_Z';                phypar{15,2} = [6    0 0 0 0];
phypar{16,1} = 'L';                         phypar{16,2} = 2.0d+1;
phypar{17,1} = 'recycling';                 phypar{17,2} = 1.00d+0;
phypar{18,1} = 'sintheta';                  phypar{18,2} = 1.0d-1;
phypar{19,1} = 'case_AMJUEL';               phypar{19,2} = '.true.';
phypar{20,1} = 'minimum_temperature';       phypar{20,2} = 0.1d-0;
phypar{21,1} = 'gas_puff_source';           phypar{21,2} = 0.0d+24;
phypar{22,1} = 'gas_puff_location';         phypar{22,2} = 5.5d+0;
phypar{23,1} = 'gas_puff_width';            phypar{23,2} = 1.0d+20;
phypar{24,1} = 'elm_start_time';            phypar{24,2} = 1000;
phypar{25,1} = 'elm_ramp_time';             phypar{25,2} = 200;
phypar{26,1} = 'elm_time_between';          phypar{26,2} = 2d+3;
phypar{27,1} = 'elm_expelled_heat';         phypar{27,2} = 40d+3;
phypar{28,1} = 'elm_expelled_particles';    phypar{28,2} = 5d+18;
phypar{29,1} = 'switch_elm_heat_flux';      phypar{29,2} = 0;
phypar{30,1} = 'switch_elm_density';        phypar{30,2} = 0;
phypar{31,1} = 'switch_elm_series';         phypar{31,2} = 0;
phypar{32,1} = 'radial_loss_factor';        phypar{32,2} = 0.0d+00;
phypar{33,1} = 'radial_loss_gaussian';      phypar{33,2} = 0;
phypar{34,1} = 'radial_loss_width';         phypar{34,2} = 1.0d+20;
phypar{35,1} = 'radial_loss_location';      phypar{35,2} = 0.0d+00;
phypar{36,1} = 'density_ramp_rate';         phypar{36,2} = 0.0d+00;
        

    case 'v4.0.1'
numpar = struct;
    
    %% Numerical parameters ($div1d_numerics)
numpar{01,1} = 'restart';                   numpar{01,2} = '.true.';
numpar{02,1} = 'Nx';                        numpar{02,2} = 500;
numpar{03,1} = 'dxmin';                     numpar{03,2} = 0.05d-0;
numpar{04,1} = 'ntime';                     numpar{04,2} = 15000;
numpar{05,1} = 'nout';                      numpar{05,2} = 100;
numpar{06,1} = 'delta_t';                   numpar{06,2} = 1.0d-6;
numpar{07,1} = 'abstol';                    numpar{07,2} = 1.0d-6;
numpar{08,1} = 'reltol';                    numpar{08,2} = 1.0d-6;
numpar{09,1} = 'method';                	numpar{09,2} = 227;
numpar{10,1} = 'max_step';                  numpar{10,2} = 100000;
numpar{11,1} = 'max_attempts';              numpar{11,2} = 1000;
numpar{12,1} = 'viscosity';                 numpar{12,2} = 4.0d-2;
numpar{13,1} = 'evolve_density';            numpar{13,2} = 1;
numpar{14,1} = 'evolve_momentum';           numpar{14,2} = 1;
numpar{15,1} = 'evolve_energy';             numpar{15,2} = 1;
numpar{16,1} = 'evolve_neutral';            numpar{16,2} = 1;
numpar{17,1} = 'switch_density_source';     numpar{17,2} = 1.0d+0;
numpar{18,1} = 'switch_momentum_source';    numpar{18,2} = 1.0d+0;
numpar{19,1} = 'switch_energy_source';      numpar{19,2} = 1.0d+0;
numpar{20,1} = 'switch_neutral_source';     numpar{20,2} = 1.0d+0;
numpar{21,1} = 'switch_charge_exchange';    numpar{21,2} = 1.0d+0;
numpar{22,1} = 'switch_ionization';         numpar{22,2} = 1.0d+0;
numpar{23,1} = 'switch_excitation';         numpar{23,2} = 1.0d+0;
numpar{24,1} = 'switch_convective_heat';    numpar{24,2} = 1.0d+0;
numpar{25,1} = 'switch_recombination';      numpar{25,2} = 1.0d+0;
numpar{26,1} = 'switch_impurity_radiation'; numpar{26,2} = 1.0d+0;
numpar{27,1} = 'simple_sol';                numpar{27,2} = '.false.';
numpar{28,1} = 'filter_sources';            numpar{28,2} = '.false.';
%% Physics parameters ($div1d_physics)
phypar{01,1} = 'initial_n';                 phypar{01,2} = 3.1d+19;
phypar{02,1} = 'initial_v';                 phypar{02,2} = 0.0d+5;
phypar{03,1} = 'initial_T';                 phypar{03,2} = 1.044d+2;
phypar{04,1} = 'initial_a';                 phypar{04,2} = 1.0d+13;
phypar{05,1} = 'gamma';                     phypar{05,2} = 6.0d+0;
phypar{06,1} = 'mass';                      phypar{06,2} = 3.3436d-27;
phypar{07,1} = 'q_parX';                    phypar{07,2} = 0.5d+8;
phypar{08,1} = 'flux_expansion';            phypar{08,2} = 1.0d+0; 
%  B_field = 1/ ( 1 + (flux_expansion-1)*x/L)
phypar{09,1} = 'Gamma_X';                   phypar{09,2} = 0.0d+24;
phypar{10,1} = 'energy_loss_ion';           phypar{10,2} = 3.0d+1;
phypar{11,1} = 'neutral_residence_time';    phypar{11,2} = 1d+20;
phypar{12,1} = 'redistributed_fraction';    phypar{12,2} = 0.0d+0;
phypar{13,1} = 'num_impurities';            phypar{13,2} = 5;
phypar{14,1} = 'impurity_concentration';    phypar{14,2} = [0.01 0 0 0 0];
phypar{15,1} = 'impurity_Z';                phypar{15,2} = [6    0 0 0 0];
phypar{16,1} = 'L';                         phypar{16,2} = 2.0d+1;
phypar{17,1} = 'recycling';                 phypar{17,2} = 1.00d+0;
phypar{18,1} = 'sintheta';                  phypar{18,2} = 1.0d-1;
phypar{19,1} = 'case_AMJUEL';               phypar{19,2} = '.true.';
phypar{20,1} = 'minimum_temperature';       phypar{20,2} = 0.1d-0;
phypar{21,1} = 'gas_puff_source';           phypar{21,2} = 0.0d+24;
phypar{22,1} = 'gas_puff_location';         phypar{22,2} = 5.5d+0;
phypar{23,1} = 'gas_puff_width';            phypar{23,2} = 1.0d+20;
phypar{24,1} = 'elm_start_time';            phypar{24,2} = 1000;
phypar{25,1} = 'elm_ramp_time';             phypar{25,2} = 200;
phypar{26,1} = 'elm_time_between';          phypar{26,2} = 2d+3;
phypar{27,1} = 'elm_expelled_heat';         phypar{27,2} = 40d+3;
phypar{28,1} = 'elm_expelled_particles';    phypar{28,2} = 5d+18;
phypar{29,1} = 'switch_elm_heat_flux';      phypar{29,2} = 0;
phypar{30,1} = 'switch_elm_density';        phypar{30,2} = 0;
phypar{31,1} = 'switch_elm_series';         phypar{31,2} = 0;
phypar{32,1} = 'radial_loss_factor';        phypar{32,2} = 0.0d+00;
phypar{33,1} = 'radial_loss_gaussian';      phypar{33,2} = 0;
phypar{34,1} = 'radial_loss_width';         phypar{34,2} = 1.0d+20;
phypar{35,1} = 'radial_loss_location';      phypar{35,2} = 0.0d+00;
phypar{36,1} = 'density_ramp_rate';         phypar{36,2} = 0.0d+00;
phypar{37,1} = 'L_core_SOL';                phypar{37,2} = 0.0d+00;
phypar{38,1} = 'X_core_SOL';                phypar{38,2} = 0.0d+00;
phypar{39,1} = 'Gamma_core';                phypar{39,2} = 0.0d+00;
phypar{40,1} = 'Q_core';                    phypar{40,2} = 1.0d+04;
       
    case 'v4.0.2'
        numpar = struct;
        phypar = struct;

    otherwise
       
        
        disp('version has no default input specified')
end