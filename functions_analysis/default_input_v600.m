function [numpar,phypar,dat] = default_input_v600(varargin)
% Generates default input for DIV1D
% variable arguments:
%    - input_div1d('radial_losses',1,'time_dependent',1) (Default option)
%   allows to generate inputs according to previous versions of DIV1D
%   that did not have radial losses or time-dependent .dat inputs

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023
D.case = '';
P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end

switch P.case
    case "TCVgeom"
         load("TCVgeom.mat");
    otherwise
           
    numpar = struct;
    numpar.Nx = int64(500);
    numpar.ntime = int64(10000);
    numpar.dxmin = dxmin;
    numpar.delta_t = 1e-6;
    numpar.abstol = 1e-6;
    numpar.reltol = 1e-6;
    numpar.viscosity = 5;
    numpar.method  = int32(227);
    numpar.density_norm = 1e19;
    numpar.temperature_norm = 1.0;
    numpar.velocity_norm = 1e4;
    % do not evolve the background and/or core
    numpar.evolve_core = 0;
    numpar.evolve_background = 0;
    
    phypar = struct;
    phypar.gamma        = 6.0;
    phypar.L            = 18;
    phypar.L_core_SOL   = 10; 
    phypar.X_core_SOL   = 0;
    phypar.Q_core       =  270*10^3; % J/s mean(iunputs.Q_core) 
    % NOTE: less than the 360 kW going into SOLPS, goes to inner leg and core losses.
    phypar.impurity_concentration = [1 0 0 0 0]*0.03;   % round(inputs.carbon_concentration,2)
    phypar.neutral_energy = 0.5;
    phypar.initial_a = 9.51e16; 
    phypar.initial_n = 1e19;
    phypar.initial_m = phypar.initial_a;
    phypar.initial_t = 100;
     
    tau_n = 5e-5; 
    tau_m = 4e-3;
    ncore = 5*1e15;
    
    phypar.neutral_residence_time    = tau_n;
    phypar.molecule_residence_time   = tau_m;
    phypar.core_sol_neutral_ex_time  = tau_n;
    phypar.core_sol_molecule_ex_time = tau_m;
    phypar.initial_core_neutral     = ncore;
    phypar.major_radius             = 0.9; 
    phypar.sintheta                 = 0.05; 
    phypar.flux_expansion           = 1;  
    phypar.recycling                = 0.4; 
    phypar.sol_width_omp            = 0.011;
    phypar.location_omp             = 0;
    phypar.alpha_core_profile_n     = 1; 
    phypar.alpha_core_profile_q     = 1; 
    phypar.trans_expansion          = 2.1;
    phypar.mol_rec = 1.0;
    
    phypar.initial_nb = [1.186189545416666e+18 1.023472520000000e+18 3.254512813703704e+17 1.023472520000000e+18 1.186189545416666e+18];
    phypar.initial_mb = [7.348254275000001e+18 6.928278916666667e+18 1.763642667777778e+18 6.928278916666667e+18 7.348254275000001e+18];
    phypar.Gamma_core = 4.564368705685312e+21;
    
    % struct with the vectors prescribing geometry
    %datinput = 1;
    dat = struct;
    dat.B_field = nan; %geometry.B_field; 
    dat.sintheta = nan; % geometry.sintheta;
    dat.major_radius = nan; %geometry.major_radius;

end

end