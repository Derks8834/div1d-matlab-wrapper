
%% define input for DIV1D simulation
numpar = struct;
numpar.Nx = int64(Nx);
numpar.nout_steps = int64(180);
numpar.nout = int64(100);
numpar.ntime = numpar.nout*numpar.nout_steps;
numpar.dxmin = dxmin;
numpar.delta_t = 1e-6;
numpar.abstol = 1e-6;
numpar.reltol = 1e-6;
numpar.viscosity = 0.5;
numpar.method  = int32(227);
numpar.density_norm = 1e19;
numpar.temperature_norm = 1.0;
numpar.velocity_norm = 1e4;
% do not evolve the background and/or core
numpar.evolve_core = 0;
numpar.evolve_background = [0 0 0 0 0];

phypar = struct;
phypar.gamma        = 6.0;
phypar.L            = 18;
phypar.L_core_SOL   = 10;
phypar.X_core_SOL   = 0;
phypar.Q_core       =  270*10^3; % J/s mean(iunputs.Q_core)
% NOTE: less than the 360 kW going into SOLPS, goes to inner leg and core losses.
phypar.impurity_concentration = [1 0 0 0 0]*0.03;   % round(inputs.carbon_concentration,2)
phypar.neutral_energy = 4;
phypar.initial_a = 1e15;
phypar.initial_n = 1e19;
phypar.initial_m = phypar.initial_a;
phypar.initial_t = 100;
phypar.core_sol_neutral_ex = 1.0d-20; % core_ionization_fraction is used
phypar.core_sol_molecule_ex = 1.0d-20; % core_ionization_fraction is used
phypar.initial_core_neutral = 1e13;% ncore;
phypar.major_radius             = 1; %-1;
phypar.sintheta                 = 0.05; %-1;
phypar.flux_expansion           = 1; %-1;
phypar.recycling                = 0.4; % this results in a bifurcation at the target
phypar.sol_width_omp            = 0.02; %lambda_omp;
phypar.location_omp             = 0;
phypar.alpha_core_profile_n       = 0.1;
phypar.alpha_core_profile_q       = 0.1;
phypar.trans_expansion          = 1; %2.1;
phypar.mol_rec = 0.5;
phypar.sigma_nb = 0.45;

% standard setting for equilibria of external chambers that we want
phypar.core_ionization_fraction = 0.6;
phypar.core_ionization_fraction_mol = phypar.core_ionization_fraction;

tau_n = 7.5e-4;   % this we should scan eventually
tau_m = 1.5e-2;
phypar.neutral_residence_time   = tau_n;
phypar.molecule_residence_time  = tau_m;

phypar.initial_nb =  [9.181991000000000   6.876696000000000   3.384090000000000   6.876696000000000   9.181991000000000]*1e17;
phypar.initial_mb =  [4.541583000000000   2.764769000000000   1.653272000000000   2.764769000000000   4.541583000000000]*1e18;
phypar.initial_ncore = 5.13e19;

phypar.Gamma_core = 4e21;
phypar.core_volume   = 1;

%% calculated with post processing: get chamber params

phypar.initial_ncore = 5.130000000E+19;
phypar.core_confinement_time = 1.282500000E-02;
phypar.core_ext_neutral_pump = [0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 1.272807847E+03];
phypar.core_fuelling = 1.000000000E+21;
phypar.extern_neutral_ex = [0.000000000E+00 0.000000000E+00 7.379585459E+03];
phypar.extern_molecule_ex = [0.000000000E+00 0.000000000E+00 5.797949912E+02];
phypar.puff_rate_neutral =[0.000000000E+00 0.000000000E+00 0.000000000E+00 3.882982027E+21 2.640364525E+21];
phypar.puff_rate_molecule =[0.000000000E+00 0.000000000E+00 7.815651775E+18 3.000000000E+20 3.000000000E+20];
phypar.pump_rate_n =[0.000000000E+00 0.000000000E+00 1.370676213E+01 0.000000000E+00 0.000000000E+00];
phypar.pump_rate_m =[0.000000000E+00 0.000000000E+00 0.000000000E+00 4.779029182E+02 5.614367392E+02];

switch check_iteration
    case 0
        numpar.evolve_molecule = 1;
        numpar.evolve_neutral_momentum = 1;
        numpar.evolve_neutral = 1;
    case  1
        % particle conservation SOL
        % all interaction fluxes are zero for now.
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 1;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %tau_n;
        phypar.molecule_residence_time  = 1e-4; %tau_m;
        phypar.Gamma_core = 1e20;%0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.1;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;
        phypar.core_ionization_fraction = 0;
        phypar.core_ionization_fraction_mol = 0;

        % this should reach a steady state
        phypar.initial_a = 1e16;
        phypar.initial_n = 1e19;
        phypar.initial_m = 1e18; %phypar.initial_a;
        % with a sum of these densities times the volume of the cells
        % (1e19 + 1e15 +2e15)*V = ??

    case 2
        % particle conservation without extern flows and without molecules
        % all internal balances are active now
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 0;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e24; %1e20; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 0.0;
        phypar.recycling = 1;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = 1.0; %phypar.initial_a;
        % with a sum of these densities times the volume of the cells
        % (1e19 + 1e15 +2e15)*V = ??
    case 3
        % particle conservation with some extern flows
        % all internal balances are active now
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 0;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %1e20; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 1e20; %0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0; % no molecule recycling
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = 1.0;% phypar.initial_a;
        % with a sum of these densities times the volume of the cells
        % (1e19 + 1e15 +2e15)*V = ??
        % seems to go ok, lets try without external sources

    case 4
        % turn on flux expansion gives me geometric error
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 0;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e24; %1e20; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 0.0; %1e20; %0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0; % no molecule recycling
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 2;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = 1.0;% phypar.initial_a;
    case 5
        % turn on flux expansion gives me geometric error but give a way to
        % balance errors
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 0;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %1e20; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 1e20; %0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0; % no molecule recycling
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 2;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = 1.0;% phypar.initial_a;
        % with a sum of these densities times the volume of the cells
        % (1e19 + 1e15 +2e15)*V = ??
        % seems to go ok, lets try without external sources

    case 6 % start looking at the molecules
        % add target molecule recycling
        % (this is the only source without external reservoirs)
        % add flux from core, because the molecuels that recycle are sink now
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 0;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %1e20; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 1e20; %0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.1; % no molecule recycling
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;


        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = 1.0;

    case 7 % start looking at the molecules
        % add target molecule recycling
        % (this is the only source without external reservoirs)
        % add flux from core, because the molecuels that recycle are sink now

        % particle conservation SOL
        % all interaction fluxes are zero for now.
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 1;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %tau_n;
        phypar.molecule_residence_time  = 1e-3; %tau_m;
        phypar.Gamma_core = 1e20;%0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.5;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;
        phypar.core_ionization_fraction = 0;
        phypar.core_ionization_fraction_mol = 0;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = phypar.initial_a;

        numpar.switch_momentum_transfer_atoms = 0.0; %d0 ! multiplier of atomic momentum transfer rate
        numpar.switch_momentum_transfer_molecules = 1.0;%d0 ! multiplier of molecular momentum transfer rate

        numpar.switch_dissociation_molecules = 0.0;%d0 ! multiplier of molecular dissociation rate
        numpar.switch_ionization_molecules = 0.0;%d0 ! multiplier of molecular ionization rate
        numpar.switch_charge_exchange_molecules = 0.0;%d0 ! multiplier of molecular charge exchange rate
        numpar.switch_dissociation_H2plus = 0.0;%d0 ! multiplier of H2plus dissociation rate
        numpar.switch_dissociative_recombination_H2plus = 0.0;%d0 ! multiplier of H2+ dissociative recombination rate
        numpar.switch_dissociative_ionization_H2plus = 0.0;%d0 ! multiplier of H2+ dissociative ionization rate
        numpar.switch_energy_molecules = 0.0;%d0 ! multiplier of molecular energy loss rate
        numpar.switch_dissociative_attachment = 0.0; %d0 ! multiplier of dissociative attachment rate
        numpar.switch_charge_exchange_Hmin = 0.0;%d0 ! multiplier of charge exchange rate on Hmin
        numpar.switch_ionization_Hmin = 0.0; %d0
        % does not run

    case 8 % start looking at the molecules
        % add target molecule recycling
        % (this is the only source without external reservoirs)
        % add flux from core, because the molecuels that recycle are sink now

        % particle conservation SOL
        % all interaction fluxes are zero for now.
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 1;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 1e20;%0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.1;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;
        phypar.core_ionization_fraction = 0;
        phypar.core_ionization_fraction_mol = 0;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = phypar.initial_a;

      
        numpar.switch_momentum_transfer_atoms = 1.0; %d0 ! multiplier of atomic momentum transfer rate
        numpar.switch_momentum_transfer_molecules = 1.0;%d0 ! multiplier of molecular momentum transfer rate

        numpar.switch_dissociation_molecules = 1.0;%d0 ! multiplier of molecular dissociation rate
        numpar.switch_ionization_molecules = 0.0;%d0 ! multiplier of molecular ionization rate
        numpar.switch_charge_exchange_molecules = 0.0;%d0 ! multiplier of molecular charge exchange rate
        numpar.switch_dissociation_H2plus = 1.0;%d0 ! multiplier of H2plus dissociation rate
        numpar.switch_dissociative_recombination_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative recombination rate
        numpar.switch_dissociative_ionization_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative ionization rate

        numpar.switch_energy_molecules = 1.0;%d0 ! multiplier of molecular energy loss rate

        numpar.switch_dissociative_attachment = 1.0; %d0 ! multiplier of dissociative attachment rate
        numpar.switch_charge_exchange_Hmin = 0.0;%d0 ! multiplier of charge exchange rate on Hmin
        numpar.switch_ionization_Hmin = 1.0; %d0

        case 9 % start looking at the molecules
        % add target molecule recycling
        % (this is the only source without external reservoirs)
        % add flux from core, because the molecuels that recycle are sink now

        % particle conservation SOL
        % all interaction fluxes are zero for now.
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 1;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 1e20;%0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.0;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;
        phypar.core_ionization_fraction = 0;
        phypar.core_ionization_fraction_mol = 0;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = phypar.initial_a;

       numpar.switch_molecule_source = 0.0;

        numpar.switch_momentum_transfer_atoms = 1.0; %d0 ! multiplier of atomic momentum transfer rate
        numpar.switch_momentum_transfer_molecules = 1.0;%d0 ! multiplier of molecular momentum transfer rate

        numpar.switch_dissociation_molecules = 1.0;%d0 ! multiplier of molecular dissociation rate
        numpar.switch_ionization_molecules = 1.0;%d0 ! multiplier of molecular ionization rate
        numpar.switch_charge_exchange_molecules = 1.0;%d0 ! multiplier of molecular charge exchange rate
        numpar.switch_dissociation_H2plus = 1.0;%d0 ! multiplier of H2plus dissociation rate
        numpar.switch_dissociative_recombination_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative recombination rate
        numpar.switch_dissociative_ionization_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative ionization rate
        
        numpar.switch_energy_molecules = 0.0;%d0 ! multiplier of molecular energy loss rate

        numpar.switch_dissociative_attachment = 1.0; %d0 ! multiplier of dissociative attachment rate
        numpar.switch_charge_exchange_Hmin = 1.0;%d0 ! multiplier of charge exchange rate on Hmin
        numpar.switch_ionization_Hmin = 1.0; %d0

           case 10 % start looking at the molecules
        % add target molecule recycling
        % (this is the only source without external reservoirs)
        % add flux from core, because the molecuels that recycle are sink now

        % particle conservation SOL
        % all interaction fluxes are zero for now.
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 1;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %tau_n;
        phypar.molecule_residence_time  = 1e-4; %tau_m;
        phypar.Gamma_core = 1e20;%0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.1;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;
        phypar.core_ionization_fraction = 0;
        phypar.core_ionization_fraction_mol = 0;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = 4*1e18; %phypar.initial_a;

        numpar.switch_molecule_source = 0.0;
        % below cannot be zero because divide by 0 in Diffusion calcs.
        numpar.switch_momentum_transfer_atoms = 1.0; %d0 ! multiplier of atomic momentum transfer rate
        numpar.switch_momentum_transfer_molecules = 1.0;%d0 ! multiplier of molecular momentum transfer rate
        
        numpar.switch_dissociation_molecules = 0.0;%d0 ! multiplier of molecular dissociation rate

        % below cannot both be zero
        numpar.switch_ionization_molecules = 1.0;%d0 ! multiplier of molecular ionization rate
        numpar.switch_charge_exchange_molecules = 0.0;%d0 ! multiplier of molecular charge exchange rate
         
        % below cannot all be zero
        numpar.switch_dissociation_H2plus = 1.0;%d0 ! multiplier of H2plus dissociation rate
        numpar.switch_dissociative_recombination_H2plus = 0.0;%d0 ! multiplier of H2+ dissociative recombination rate
        numpar.switch_dissociative_ionization_H2plus = 0.0;%d0 ! multiplier of H2+ dissociative ionization rate
        
        numpar.switch_energy_molecules = 0.0;%d0 ! multiplier of molecular energy loss rate

        numpar.switch_dissociative_attachment = 1.0; % turn off Hmin %d0 ! multiplier of dissociative attachment rate
        % below cannot both be zero
        numpar.switch_charge_exchange_Hmin = 1.0;%d0 ! multiplier of charge exchange rate on Hmin
        numpar.switch_ionization_Hmin = 0.0; %d0

           case 11 % start looking at the molecules
        % add target molecule recycling
        % (this is the only source without external reservoirs)
        % add flux from core, because the molecuels that recycle are sink now

        % particle conservation SOL
        % all interaction fluxes are zero for now.
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 1;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 1e20;%0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.1;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;
        phypar.core_ionization_fraction = 0;
        phypar.core_ionization_fraction_mol = 0;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = phypar.initial_a;

        numpar.switch_molecule_source = 0.0;
        % below cannot be zero because divide by 0 in Diffusion calcs.
        numpar.switch_momentum_transfer_atoms = 1.0; %d0 ! multiplier of atomic momentum transfer rate
        numpar.switch_momentum_transfer_molecules = 1.0;%d0 ! multiplier of molecular momentum transfer rate
        
        numpar.switch_dissociation_molecules = 0.0;%d0 ! multiplier of molecular dissociation rate

        % below cannot both be zero
        numpar.switch_ionization_molecules = 1.0;%d0 ! multiplier of molecular ionization rate
        numpar.switch_charge_exchange_molecules = 0.0;%d0 ! multiplier of molecular charge exchange rate
         
        % below cannot all be zero
        numpar.switch_dissociation_H2plus = 1.0;%d0 ! multiplier of H2plus dissociation rate
        numpar.switch_dissociative_recombination_H2plus = 0.0;%d0 ! multiplier of H2+ dissociative recombination rate
        numpar.switch_dissociative_ionization_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative ionization rate
        
        numpar.switch_energy_molecules = 0.0;%d0 ! multiplier of molecular energy loss rate

        numpar.switch_dissociative_attachment = 1.0; % turn off Hmin %d0 ! multiplier of dissociative attachment rate
        % below cannot both be zero
        numpar.switch_charge_exchange_Hmin = 0.0;%d0 ! multiplier of charge exchange rate on Hmin
        numpar.switch_ionization_Hmin = 1.0; %d0

           case 12 % start looking at the molecules
        % add target molecule recycling
        % (this is the only source without external reservoirs)
        % add flux from core, because the molecuels that recycle are sink now

        % particle conservation SOL
        % all interaction fluxes are zero for now.
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 1;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %tau_n;
        phypar.molecule_residence_time  = 1e-4; %tau_m;
        phypar.Gamma_core = 1e20;%0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.1;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;
        phypar.core_ionization_fraction = 0;
        phypar.core_ionization_fraction_mol = 0;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = phypar.initial_a;

        numpar.switch_molecule_source = 0.0;
        % below cannot be zero because divide by 0 in Diffusion calcs.
        numpar.switch_momentum_transfer_atoms = 1.0; %d0 ! multiplier of atomic momentum transfer rate
        numpar.switch_momentum_transfer_molecules = 1.0;%d0 ! multiplier of molecular momentum transfer rate
        
        numpar.switch_dissociation_molecules = 0.0;%d0 ! multiplier of molecular dissociation rate

        % below cannot both be zero
        numpar.switch_ionization_molecules = 1.0;%d0 ! multiplier of molecular ionization rate
        numpar.switch_charge_exchange_molecules = 0.0;%d0 ! multiplier of molecular charge exchange rate
         
        % below cannot all be zero
        numpar.switch_dissociation_H2plus = 0.0;%d0 ! multiplier of H2plus dissociation rate
        numpar.switch_dissociative_recombination_H2plus = 0.0;%d0 ! multiplier of H2+ dissociative recombination rate
        numpar.switch_dissociative_ionization_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative ionization rate
        
        numpar.switch_energy_molecules = 0.0;%d0 ! multiplier of molecular energy loss rate

        numpar.switch_dissociative_attachment = 0.0; % turn off Hmin %d0 ! multiplier of dissociative attachment rate
        % below cannot both be zero
        numpar.switch_charge_exchange_Hmin = 1.0;%d0 ! multiplier of charge exchange rate on Hmin
        numpar.switch_ionization_Hmin = 0.0; %d0

           case 13 % start looking at the molecules
        % add target molecule recycling
        % (this is the only source without external reservoirs)
        % add flux from core, because the molecuels that recycle are sink now

        % particle conservation SOL
        % all interaction fluxes are zero for now.
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 1;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 1e20;%0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.1;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 1;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;
        phypar.core_ionization_fraction = 0;
        phypar.core_ionization_fraction_mol = 0;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = phypar.initial_a;

        numpar.switch_molecule_source = 1.0;
        % below cannot be zero because divide by 0 in Diffusion calcs.
        numpar.switch_momentum_transfer_atoms = 1.0; %d0 ! multiplier of atomic momentum transfer rate
        numpar.switch_momentum_transfer_molecules = 1.0;%d0 ! multiplier of molecular momentum transfer rate
        
        numpar.switch_dissociation_molecules = 1.0;%d0 ! multiplier of molecular dissociation rate

        % below cannot both be zero
        numpar.switch_ionization_molecules = 1.0;%d0 ! multiplier of molecular ionization rate
        numpar.switch_charge_exchange_molecules = 1.0;%d0 ! multiplier of molecular charge exchange rate
         
        % below cannot all be zero
        numpar.switch_dissociation_H2plus = 1.0;%d0 ! multiplier of H2plus dissociation rate
        numpar.switch_dissociative_recombination_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative recombination rate
        numpar.switch_dissociative_ionization_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative ionization rate
        
        numpar.switch_energy_molecules = 1.0;%d0 ! multiplier of molecular energy loss rate

        numpar.switch_dissociative_attachment = 1.0; % turn off Hmin %d0 ! multiplier of dissociative attachment rate
        % below cannot both be zero
        numpar.switch_charge_exchange_Hmin = 1.0;%d0 ! multiplier of charge exchange rate on Hmin
        numpar.switch_ionization_Hmin = 1.0; %d0
    case 14 % B field and molecules
         % start looking at the molecules
        % add target molecule recycling
        % (this is the only source without external reservoirs)
        % add flux from core, because the molecuels that recycle are sink now

        % particle conservation SOL
        % all interaction fluxes are zero for now.
        numpar.evolve_core = 0;
        numpar.evolve_background = [0 0 0 0 0];
        numpar.evolve_molecule = 1;
        % turn off external influxes of particles
        phypar.neutral_residence_time   = 1e-4; %tau_n;
        phypar.molecule_residence_time  = 1e24; %tau_m;
        phypar.Gamma_core = 1e20;%0.0;
        phypar.recycling = 1;
        phypar.mol_rec = 0.1;
        phypar.trans_expansion = 1;
        phypar.flux_expansion  = 2;
        phypar.major_radius = 1;
        phypar.sintheta = 0.05;
        phypar.core_ionization_fraction = 0;
        phypar.core_ionization_fraction_mol = 0;

        % this should reach a steady state
        phypar.initial_a = 1e15;
        phypar.initial_n = 1e19;
        phypar.initial_m = phypar.initial_a;

        numpar.switch_molecule_source = 1.0;
        % below cannot be zero because divide by 0 in Diffusion calcs.
        numpar.switch_momentum_transfer_atoms = 1.0; %d0 ! multiplier of atomic momentum transfer rate
        numpar.switch_momentum_transfer_molecules = 1.0;%d0 ! multiplier of molecular momentum transfer rate
        
        numpar.switch_dissociation_molecules = 1.0;%d0 ! multiplier of molecular dissociation rate

        % below cannot both be zero
        numpar.switch_ionization_molecules = 1.0;%d0 ! multiplier of molecular ionization rate
        numpar.switch_charge_exchange_molecules = 1.0;%d0 ! multiplier of molecular charge exchange rate
         
        % below cannot all be zero
        numpar.switch_dissociation_H2plus = 1.0;%d0 ! multiplier of H2plus dissociation rate
        numpar.switch_dissociative_recombination_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative recombination rate
        numpar.switch_dissociative_ionization_H2plus = 1.0;%d0 ! multiplier of H2+ dissociative ionization rate
        
        numpar.switch_energy_molecules = 1.0;%d0 ! multiplier of molecular energy loss rate

        numpar.switch_dissociative_attachment = 1.0; % turn off Hmin %d0 ! multiplier of dissociative attachment rate
        % below cannot both be zero
        numpar.switch_charge_exchange_Hmin = 1.0;%d0 ! multiplier of charge exchange rate on Hmin
        numpar.switch_ionization_Hmin = 1.0; %d0
end

% calculate target angle
dr = abs(0.887 - 0.868);
dz = abs(-0.75 +0.65);
ang = rad2deg(atan(dz/dr)); % gives 79 deg (we use 80)
phypar.pol_target_angle =  [ 90.0  100.0 ]; % 0 to 1,4, 180 to 2 and 5

% struct with the vectors prescribing geometry
datinput = 1;
dat.B_field = geometry.B_field;
dat.sintheta = geometry.sintheta;
dat.major_radius = geometry.major_radius;
dat.major_height = geometry.major_height;
dat.romrz = geometry.romrz;
dat.domain_bound = geometry.domain_bound;
dat.vesrz = geometry.vesrz;
dat.ir = geometry.ir;
dat.or = geometry.or;
dat.iz = geometry.iz;
dat.oz = geometry.oz;
dat.sol_normal = geometry.normal_vector;