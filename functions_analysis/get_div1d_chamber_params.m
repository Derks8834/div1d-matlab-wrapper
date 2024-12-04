function [sim,set] = get_div1d_chamber_params(in,out,varargin)
%GET_DIV1D_CHAMBER_PARAMS 
% example: [sim,set] = get_div1d_chamber_params(in,out, pump_n, pump_m, tcore,ncore)
%
%  inputs: follow from div1d_read.m
%   - in = div1d input struct
%   - out = div1d1 output struct
% D = struct;
% core losses
% D.core_density = 0;
% D.tau_core = 0;
% % core sources
% D.core_fuelling = 0;
% D.core_ext_molecule_pump = [0 0 0 0 0]; 
% D.core_ext_neutral_pump = [0 0 0 0 0];
% % reservoir circulating flows
% D.circflow_n_temp = 0.5; % eV
% D.leak_gaps = [0.4 0.1 0.2]; % m 
% D.leak_ato_mult = [0 0 3]; %[1 1 1];
% D.leak_mol_mult = [0 0 1/3]; %[1 1 1];
% % reservoir pumps
% 
% 
% outputs: 
%   - sim = analysis of the simulation output and the balances
%   - set = settings for the neutral chambers and core for steady state solution



%
% phypar.initial_ncore = 5.130000000E+19;
% phypar.core_confinement_time = 1.282500000E-02;
% phypar.core_ext_neutral_pump = [0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 1.247384013E+03];
% phypar.core_fuelling = 1.000000000E+21;
% phypar.extern_neutral_ex = [0 0 7.379585459E+03];
% phypar.extern_molecule_ex = [0 0 5.797949912E+02];
% phypar.puff_rate_neutral =[0.000000000E+00 0.000000000E+00 0.000000000E+00 5.471164732E+19 0.000000000E+00];
% phypar.puff_rate_molecule =[0.000000000E+00 0.000000000E+00 8.421824775E+18 3.039127229E+21 3.000000000E+20];
% phypar.pump_rate_n =[0.000000000E+00 0.000000000E+00 4.509157452E+00 0.000000000E+00 2.181648004E+03];
% phypar.pump_rate_m =[0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 6.813638781E+02];
%
%
% author g.l. derks 2024 g.l.derks@differ.nl

%constants
e_charge = 1.61e-19;

sim = struct;
set = struct;

D = struct;
% core losses
D.core_density = 0;
D.tau_core = 0;
% core sources
D.core_fuelling = 0;
D.core_ext_molecule_pump = [0 0 0 0 0]; 
D.core_ext_neutral_pump = [0 0 0 0 0];
% reservoir circulating flows
D.circflow_n_temp = 0.5; % eV
D.leak_gaps = [0.4 0.1 0.2]; % m 
D.leak_ato_mult = [0 0 3]; %[1 1 1];
D.leak_mol_mult = [0 0 1/3]; %[1 1 1];
% reservoir pumps
D.atom_puff = [0 0 0 0 0 ];
D.molecule_puff = [0 0 0 0 0];
D.atom_pump = [0 0 0 0 0];
D.mol_pump = [0 0 0 0 0];

D.core2main = 0;

% print or not
D.print = 0;

P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end

%% core balance
% dndt = -Gamma_core2sol ( fixed by given sol solution )
% + sol2core_flux + 2* sol2core_mol ( fixed by given sol solution )
% + sum_extern2core_flux + 2* sum_extern2core_mol (set core_ext_molecule/neutral_pump )
% + input.core_fuelling (set core_fuelling )

% Gamma_core2sol = (1/tau)*core_density*core_volume (set core_density or tau to match Gamma)
if P.core_density == 0  && P.tau_core > 0 % use tau-core to set core density
set.tau_core = P.tau_core;
set.core_density = set.tau_core.*in.physics.gamma_core/in.physics.core_volume;
%error('specifying tau_core not supported yet')
elseif P.core_density >0 && P.tau_core == 0 % use core density to set tau_core
set.core_density = P.core_density;
set.tau_core = set.core_density*in.physics.core_volume./in.physics.gamma_core;
else
    error(' specify either core_density or tau_core')
end


if min(P.core_ext_neutral_pump) < 0
    error('cannot specify negative neutral pump')
end
if min(P.core_ext_molecule_pump) < 0
    error('cannot specify negative molecule pump')
end

% the influx to the core should match the outflux:
if P.core_fuelling == 0 && max(P.core_ext_neutral_pump + P.core_ext_molecule_pump) > 0   % 
    disp('set core-fuelling')
    set.core_ext_molecule_pump = P.core_ext_molecule_pump;
    set.core_ext_neutral_pump = P.core_ext_neutral_pump;
    extern2core_atoms = set.core_ext_neutral_pump.*out.extern_neutral_density(end,:).* in.physics.extern_neutral_volumes;
    extern2core_molecule=set.core_ext_molecule_pump.*out.extern_molecule_density(end,:).* in.physics.extern_neutral_volumes;
    set.core_fuelling = in.physics.gamma_core ...
        -  out.sol2core_flux(end,1) - 2* out.sol2core_mol(end,1) ...
        - sum(extern2core_atoms) - 2*sum(extern2core_molecule);

elseif P.core_fuelling ~= 0 && sum(P.core_ext_neutral_pump + P.core_ext_molecule_pump) == 0 % 
     disp('set extern 2 core neutral pump')
    set.core_fuelling = P.core_fuelling;
    set.core_ext_molecule_pump = [0 0 0 0 0];
    set.core_ext_neutral_pump = [0 0 0 0 0];
    required_flux2core = in.physics.gamma_core- set.core_fuelling -  out.sol2core_flux(end,1) - 2* out.sol2core_mol(end,1);
    set.core_ext_neutral_pump(5) = required_flux2core / out.extern_neutral_density(end,5)/ in.physics.extern_neutral_volumes(5);
    % for now only room 5 provides only atoms to the core

    extern2core_atoms = set.core_ext_neutral_pump.*out.extern_neutral_density(end,:).*in.physics.extern_neutral_volumes;
    extern2core_molecule=set.core_ext_molecule_pump.*out.extern_molecule_density(end,:).* in.physics.extern_neutral_volumes;
else
    error(' specify either core_fuelling or core_ext_neutral_pump')
end

%% reservoir balance
% for both atoms and molecules
% dndt =     sol2extern +  tar2extern % ( fixed by given sol solution )
       %  -  extern2core              % ( fixed by above core closure )
       %  +  sol2extern_ion           % ( fixed by given sol solution )
       %  +  extern flows             % ( set by extern_atom/molec_ex )
       %  +  puff                     % ( set by puff_rate_atom/molec )
       %  +  pump                     % ( set by pump_rate_n/m        )

% circulating flows
vato = (e_charge*P.circflow_n_temp *8./ in.physics.mass./pi)^0.5;
vmol = (e_charge*P.circflow_n_temp *8./ in.physics.mass./2./pi)^0.5;

if (max(P.leak_ato_mult) >= 0) % thermal 1-way maxwellian 1/4 n*v *A
set.extern_neutral_ex = 2*pi*P.leak_gaps*1/4*vato.*P.leak_ato_mult;
set.extern_molecule_ex = 2*pi*P.leak_gaps*1/4*vmol.*P.leak_mol_mult;
elseif (min(P.leak_ato_mult < 0))
    error('leakage factors are input and not calculated in this function')
end

%	! fluxes go as : (1) from 5to1, (2) from 2to3, (3) from 3to4
% extern_fluxes =          (/(extern_neutral(5)-extern_neutral(1))*extern_ex(1), &   ! from 5->1
%                          (extern_neutral(2)-extern_neutral(3))*extern_ex(2), &	    ! from 2->3
%                          (extern_neutral(3)-extern_neutral(4))*extern_ex(3) /)     ! from 3->4
% then sources are calculated as:
% ext_src = (/ extern_neutral_flux(1) , & ! positive direction for neutral flows outside the plasma is anti-clockwise
% 			  -extern_neutral_flux(2) , &
% 			   extern_neutral_flux(2) - extern_neutral_flux(3) , &
% 			   extern_neutral_flux(3) , &
% 			  -extern_neutral_flux(1) /)	

% atoms 
ext_ato_flw = [set.extern_neutral_ex(1)*( out.extern_neutral_density(end,5)- out.extern_neutral_density(end,1))  ...
      set.extern_neutral_ex(2)*( out.extern_neutral_density(end,2)- out.extern_neutral_density(end,3)) ...
      set.extern_neutral_ex(3)*( out.extern_neutral_density(end,3)- out.extern_neutral_density(end,4)) ];
ext_ato_flw_src(1:5) =[ext_ato_flw(1), ...
                            -ext_ato_flw(2),  ...
                            ext_ato_flw(2)-ext_ato_flw(3), ...
                            ext_ato_flw(3), ...
                            -ext_ato_flw(1)]; % D/s
% molecules    
ext_mol_flw = [set.extern_molecule_ex(1)*( out.extern_molecule_density(end,5)- out.extern_molecule_density(end,1))... 
      set.extern_molecule_ex(2)*( out.extern_molecule_density(end,2)- out.extern_molecule_density(end,3)) ...
      set.extern_molecule_ex(3)*( out.extern_molecule_density(end,3)- out.extern_molecule_density(end,4)) ];
ext_mol_flw_src(1:5)=[ext_mol_flw(1), ...
                            -ext_mol_flw(2),  ...
                            ext_mol_flw(2)-ext_mol_flw(3), ...
                            ext_mol_flw(3), ...
                            -ext_mol_flw(1)]; % D2/s

%consider situation without any puffing!
set.atom_puff = P.atom_puff; % *inputs.puff(solps_seilect); % D/s
set.molecule_puff =  P.molecule_puff; %*inputs.puff(solps_select); %D2/s


% (residuals)
set.ext_res_ato(1:5) =               out.sol2extern_flux(end,:) ... 
                                    +  out.tar2extern_flux(end,:) ...
                                    + [0 0  out.sum_sol2extern_ion_flux(end)+P.core2main 0 0]...
                                    + ext_ato_flw_src ...
                                    - extern2core_atoms ...                                        
                                    + set.atom_puff; % no neutrals puffed in SOLPS;
set.ext_res_mol(1:5) =               out.sol2extern_mol(end,:) ... 
                                    +  out.tar2extern_mol(end,:) ...
                                    -  out.extern2core_mol(end,:) ...
                                    + [0 0  out.sum_sol2extern_ion_mol(end) 0 0]...
                                    + ext_mol_flw_src ... 
                                    - extern2core_molecule ...
                                    + set.molecule_puff; % the puffing in SOLPS should be used here;

% check total balance: this includes everything
set.totalresidual = sum(set.ext_res_ato(1:5) + 2.0*set.ext_res_mol(1:5));
            
for i_room = 1:5
if set.ext_res_ato(i_room) < 0
    set.puff_ato_closure(i_room) = set.ext_res_ato(i_room)*-1;
    set.pump_rate_n(i_room) = 0.0;
else
    set.puff_ato_closure(i_room) = 0.0;
    set.pump_rate_n(i_room) = set.ext_res_ato(i_room)./ out.extern_neutral_density(end,i_room)./ in.physics.extern_neutral_volumes(i_room);
end
if set.ext_res_mol(i_room) < 0
    set.puff_mol_closure(i_room) = set.ext_res_mol(i_room)*-1;
    set.pump_rate_m(i_room) = 0.0;
else
    set.puff_mol_closure(i_room) = 0.0;
    set.pump_rate_m(i_room) = set.ext_res_mol(i_room)./ out.extern_molecule_density(end,i_room)./ in.physics.extern_neutral_volumes(i_room);
end


end

if P.print == 1
% settings to use for core:
fprintf('phypar.initial_ncore = %4.9E;\n',set.core_density);
fprintf('phypar.core_confinement_time = %4.9E;\n',set.tau_core);
fprintf('phypar.core_ext_neutral_pump = [%4.9E %4.9E %4.9E %4.9E %4.9E];\n', set.core_ext_neutral_pump);
fprintf('phypar.core_fuelling = %4.9E;\n',set.core_fuelling);
% settings to use for reservoirs:
fprintf('phypar.extern_neutral_ex = [%4.9E %4.9E %4.9E];\n', set.extern_neutral_ex(:));
fprintf('phypar.extern_molecule_ex = [%4.9E %4.9E %4.9E];\n', set.extern_molecule_ex(:));
tmp = set.puff_ato_closure+set.atom_puff;
fprintf('phypar.puff_rate_neutral =[%4.9E %4.9E %4.9E %4.9E %4.9E];\n', tmp(:) );
% of which X is input
fprintf('%phypar.puff_rate_neutral_prescribed =[%4.9E %4.9E %4.9E %4.9E %4.9E];\n', set.atom_puff(:) );
tmp = set.puff_mol_closure+set.molecule_puff;
fprintf('phypar.puff_rate_molecule =[%4.9E %4.9E %4.9E %4.9E %4.9E];\n', tmp(:) );
% of which X is input
fprintf('%phypar.puff_rate_molecule_prescribed =[%4.9E %4.9E %4.9E %4.9E %4.9E];\n', set.molecule_puff(:) );
fprintf('phypar.pump_rate_n =[%4.9E %4.9E %4.9E %4.9E %4.9E];\n', set.pump_rate_n(:));
fprintf('phypar.pump_rate_m =[%4.9E %4.9E %4.9E %4.9E %4.9E];\n', set.pump_rate_m(:));   
end

end

