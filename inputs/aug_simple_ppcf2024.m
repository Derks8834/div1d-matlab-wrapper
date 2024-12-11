numpar = struct;
numpar.Nx = int64(500);
numpar.ntime = int64(20000);
numpar.dxmin = dxmin;
numpar.delta_t = 1e-6;
numpar.abstol = 1e-6;
numpar.reltol = 1e-6;
numpar.viscosity = 5;
numpar.method  = int32(227);
numpar.density_norm = 1e19;
numpar.temperature_norm = 1.0;
numpar.velocity_norm = 1e4;
numpar.evolve_core = 0;
numpar.evolve_background = [0 0 0 0 0];

phypar = struct;
phypar.alpha_core_profile_n = 1; 
phypar.gamma        = 6.0;
phypar.L            = 24.1;
phypar.L_core_SOL   = L-9.3; 
phypar.X_core_SOL   = 0;
phypar.flux_expansion = 1.08; 
phypar.impurity_concentration = [1 0 0 0 0]*0.00; 
phypar.neutral_energy = 5; % the 0.5 eV in Derks PPCF2024 is too low.

phypar.redistributed_fraction = 0; 
phypar.switch_imp_distribution = int64(0);
phypar.initial_a = mean(inputs.neutral_density);
phypar.initial_n = 1e19;
phypar.initial_t = 100;


phypar.initial_nb = [1 1 6.74 6.57 4.62e2]*1e17;
phypar.initial_mb = phypar.initial_nb;
phypar.Gamma_core   = 0.9e22; 
phypar.Q_core       = 750e3;

phypar.neutral_residence_time   = 1/555; 
phypar.core_ionization_fraction = 0.88;
phypar.major_radius             = 1.8; 
phypar.recycling                = 0.5;
phypar.sol_width_omp            = 0.019; 
phypar.sintheta                 = 0.1;
phypar.alpha_core_profile_q       = 1;
phypar.trans_expansion          = 1.01;
phypar.location_omp = 0;


phypar.initial_ncore = 4.780000000E+19;
phypar.core_confinement_time = 6.736841689E-02;
phypar.core_ext_neutral_pump = [0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 1.00E+02];
phypar.core_fuelling = 1.000000000E+01;
phypar.extern_neutral_ex = [0.000000000E+00 0.000000000E+00 7.379585459E+03];
phypar.extern_molecule_ex = [0.000000000E+00 0.000000000E+00 5.797949912E+02];
phypar.puff_rate_neutral =[0.000000000E+00 0.000000000E+00 0.000000000E+00 1.0E+20 1.0E+20];
phypar.puff_rate_molecule =[0.000000000E+00 0.000000000E+00 0.000000000E+00 1.0E+20 1.0E+20];
phypar.pump_rate_n =[0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00];
phypar.pump_rate_m =[0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00];
   
    
