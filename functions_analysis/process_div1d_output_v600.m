function [out] = process_div1d_output_v600(o,in,varargin)
% [tmp]= process_div1d_output_1_dev(div1doutput,div1dinput,'solps',solps,'plot',1,'Ntime',timestamp,...
%      'printerror',0,'print',1,'fignum',ishot*10,'title',string,...
%      'Dtitle',string,'avg',1,'cratio',P.cratio,...);
%      'save',saveit,'allpossiblefigures',1);
%
%Tests on DIV1D solution, original from Egbert Westerhof 2016-2020
%  1. Perform checks on stationarity and accuracy of the solution
%  2. Calculate 2PM equivalent loss factors for momentum and energy from DIV1D
%  3. comparison with SOLPS 1D interpreted profile possible
%
% INPUT:
%   - div1doutput, div1dinput, taken from resultsfolder with extracted div1d data
%
% PARAMETERS
%   'print'         ==0 no print ==1 print comparison
%   'fignum'        ==16 figure number for the overview errors
%   'Gamma'         ==6 default, sheath heat transmission (not standard in div1d output)
%   'LineWidth'     ==1 default
%   'FontSize'      ==13 default
%   'dxmin'         ==0.05 default minimum grid size\
%   'plot'          ==0 default, ==1 plot figures
%   'hold'          ==0 default, ==1 hold on;
%   'solps'         ==0 default, == struct with arrays equivalent arrays to div1d
%
% OUTPUT:
%   out = struct;
% contains processed outputs
%   out.pressure = [DIV1D];                  % neutral flux [0 L]
%   out.L10eV    = L10eV; 
%   out.S2PM = S2PM;
%   out.D2PM = D2PM;
%   out.Derror = error; == error on balance equations
%   out.SDerror = 1;
%   X2PM.f_mom; X2PM.f_pow  == 2 point model calculations for SOLPS-DIV1D
%   X2PM.T_u = [2PMfit, 2PMK-R,   DIV1D];    % Temperatures upstream
%   X2PM.T_t = [2PMfit, 2PMK-R,   DIV1D];    % Temperatures target
%   X2PM.n_t = [2PMfit, 2PMK-R,   DIV1D];    % densities target

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

%% Handling parameters
% Default paramters
% -- physics ---
D.gamma = in.physics.gamma;
D.cratio = 0.03;
% -- numerics ---
%D.gamma = input.numerics.6;
D.dxmin = 0.05; 
D.xlimits = [0  7.8];
D.Ntime = length(o.time);
D.solps = 0;
D.input = 0;
D.version = 'v6.0.0';

% -- plotting ---
D.plotbal = 0;
D.ploterror = 0;
D.plotinput = 0;
D.plotprofiles = 0;
D.generalposition = [100 100 200 250];
D.print = 0;
D.printerror =1;
D.plot = 0;
D.figtight = 0;
D.plot2PM = 1;
D.allpossiblefigures = 0;
D.seperatefigs =0;
D.colorindex = 0;
D.fignum = 16;
D.LineWidth = 1;
D.FontSize = 11;
D.flip = 1;
D.title = '' ;
D.Dtitle = '' ;
D.save = 0;
D.cratio = 0; 
D.hold = 0;
D.fR = 0;
D.R = 1.0;
D.avg = 0;
D.FE = 1.0;
D.tn = 1;
P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end

%% Constants
e_charge = 1.602*1e-19;
m = 1.66*10^(-27)*2;

%% initialize output
out = struct;

%% Load DIV1D Parameters
Ntime       = P.Ntime;
Nx          = length(o.X);
x           = o.X;
xcb         = o.Xcb;
L           = in.physics.l;
L_core_SOL  = in.physics.l_core_sol;
x_core_SOL  = in.physics.x_core_sol;

B_field = in.grid.b_field./in.grid.b_trans;
B_field_cb = in.grid.b_field_cb./in.grid.b_trans_cb;
B_tot = in.grid.b_field; %.*in.grid.b_trans;
B_tot_cb = in.grid.b_field_cb; %.*in.grid.b_trans_cb;

P.tn            = in.physics.neutral_residence_time;
P.tm            = in.physics.molecule_residence_time;

flux_expansion  = in.physics.flux_expansion;
if P.FE ==1
    FE2PM = flux_expansion;
    P.FE = in.physics.flux_expansion;
else
    FE2PM = P.FE; % flux expansion factor used in Two-Point Formatting
end

% calculate spacing of the grid
delta_x = diff(x); 
delta_xcb = diff(xcb);

% extract important profiles
density_cb  = interp1(x,o.density(Ntime,:),xcb,"linear", "extrap" );
temperature_cb  = interp1(x,o.temperature(Ntime,:),xcb,"linear", "extrap" );
velocity_cb     = interp1(x,o.velocity(Ntime,:),xcb,"linear", "extrap" ); 
neutral_density_cb = interp1(x,o.neutral_density(Ntime,:),xcb,"linear", "extrap" );

% momentum
total_pressure = 2.0*o.density.*o.temperature*e_charge + (o.Gamma_mom(:,2:end)+o.Gamma_mom(:,1:end-1))/2;
pressure_t     = 2.0*o.density.*o.temperature*e_charge;
pressure_ss    = pressure_t(Ntime,:);
P_kin = 2.0*o.density(Ntime,:).*o.temperature(Ntime,:)*e_charge;
P_ram = (o.Gamma_mom(Ntime,2:end)+o.Gamma_mom(Ntime,1:end-1))/2; %.*B_field
P_tot = P_kin + P_ram;

% energy
Source_Qtot = o.Source_Q(Ntime,:);
Source_Qtot(2:Nx) = Source_Qtot(2:Nx) + o.velocity(Ntime,2:Nx).*(pressure_ss(2:Nx)-pressure_ss(1:Nx-1))./delta_x(1:Nx-1);
Source_Qtot = Source_Qtot./B_tot; 

%% post process rates
for i = length(o.X):-1:1
    nenn = o.density(Ntime,i).*o.neutral_density(Ntime,i);
    rate_exc(i)     = nenn*AMJUEL_H102_R215_e_loss_excitation(o.density(Ntime,i),o.temperature(Ntime,i)).*e_charge; 
    rate_rec(i)     = o.density(Ntime,i)*o.density(Ntime,i)*AMJUEL_H46_R218_recombination(o.density(Ntime,i),o.temperature(Ntime,i));
    rate_ion_vel(i) = nenn*AMJUEL_H43_R215_ionization(o.density(Ntime,i),o.temperature(Ntime,i)).*0.5*in.physics.mass*o.velocity(Ntime,i)^2;
    rate_ion_nen(i) = nenn*AMJUEL_H43_R215_ionization(o.density(Ntime,i),o.temperature(Ntime,i)).*e_charge.*abs(in.physics.neutral_energy).*o.temperature(Ntime,i);
    rate_ion_den(i) = nenn*AMJUEL_H43_R215_ionization(o.density(Ntime,i),o.temperature(Ntime,i));
    rate_cx(i)      = nenn*AMJUEL_H219_R318_charge_exchange(o.temperature(Ntime,i),in.physics.mass); % charge exchange requires mass to scale the rate.
    rate_rec_ene(i)    = o.density(Ntime,i)*o.density(Ntime,i)*AMJUEL_104_218_e_loss_recombination(o.density(Ntime,i), o.temperature(Ntime,i));
    rate_cooling_post_c = o.density(Ntime,i)*o.density(Ntime,i)*in.physics.impurity_concentration(1)*POST_cooling_carbon(o.temperature(Ntime,i));
end
    proc.rate_cx        = rate_cx;
    proc.rate_exc       = rate_exc;
    proc.rate_rec       = rate_rec;
    proc.rate_ion_vel   = rate_ion_vel;
    proc.rate_ion_nen   = rate_ion_nen;
    proc.rate_ion_den   = rate_ion_den;
    proc.rate_rec_ene   = rate_rec_ene;
    proc.rate_cooling_post_c = rate_cooling_post_c;
    % add rates from molecules

%% Perform checks on stationarity and accuracy (DIV1D)
den = struct;
% 1. density equation
den.dgamdx = -B_tot.*diff(o.Gamma_n(Ntime,:)./B_tot_cb)./delta_xcb;
den.src = o.Source_n(Ntime,:);
den.dt = den.dgamdx + den.src;
    error = struct;
    error.bal.density = den.dt./o.density(Ntime,:);
    error.src.density = den.src - rate_ion_den + rate_rec;

mom = struct;
% 2. momentum equation
mom.bdgamdx = -B_tot.*diff(o.Gamma_mom(Ntime,:)./B_tot_cb)./delta_xcb;
pressure_cb = interp1(x,pressure_ss,xcb,'linear','extrap');
mom.dpdx = -diff(pressure_cb)./delta_xcb;
mom.visc = zeros(1,Nx);
mom.visc(2:Nx-1) = in.numerics.viscosity*(o.velocity(Ntime,3:end) + o.velocity(Ntime,1:end-2) - 2*o.velocity(Ntime,2:end-1)); 
mom.src = o.Source_v(Ntime,:);
mom.dt = mom.dpdx + mom.bdgamdx + mom.src + mom.visc;
    error.bal.momentum = mom.dt./pressure_ss; 
    error.src.momentum = mom.src + in.physics.mass*o.velocity(Ntime,1:end).*(rate_cx - rate_rec) ;

ene = struct;
% 3. energy equation
ene.bdqdx   = -B_tot.*diff(o.q_parallel(Ntime,:)./B_tot_cb)./delta_xcb;
velocity_cc = interp1(x,o.velocity(Ntime,:),x);
ene.vdpdx   = velocity_cc .* diff(pressure_cb)./delta_xcb;
ene.src     = o.Source_Q(Ntime,:);
ene.dt    = ene.bdqdx + ene.vdpdx +ene.src;
    error.bal.energy = ene.dt./(3.0.*o.density(Ntime,:).*e_charge.*o.temperature(Ntime,:));
    error.src.energy = ene.src...
                    - ( -1.5*e_charge.*o.temperature(Ntime,:).*2.*rate_rec  ...
                     -(1.5+in.physics.neutral_energy)*e_charge.*o.temperature(Ntime,:).*rate_cx ...
                     +rate_ion_nen ...
                     +rate_ion_vel ...
                     -rate_exc ...
                     +(rate_rec*13.6 - rate_rec_ene).*e_charge ...
                     - rate_cooling_post_c+ ...
                     + in.grid.e_core_source_profile_q*in.physics.q_core./in.grid.volumes);
        % this should be updated to include the molecules.

% 4. neutral equation
ato = struct;
ato.dgamadx = -diff(o.Gamma_neutral(Ntime,:))./delta_xcb;
ato.src =  o.Source_neutral(Ntime,:);
ato.dt = ato.dgamadx + ato.src;
error.bal.neutral = ato.dt./o.neutral_density(Ntime,:);
error.src.neutral = o.Source_neutral(Ntime,:)*nan; % should be completed

% 5. neutral momentum
anv = struct;
anv.dgamdx = -diff(o.Gamma_mom_neutral(Ntime,:))./delta_xcb;
pressure_nn = e_charge * neutral_density_cb * in.physics.neutral_energy; % might want to add a temperature profile to the DIV1D output ( now constant
pressure_nn_cc = e_charge * o.neutral_density(Ntime,:) * in.physics.neutral_energy;
anv.dpdx = -diff(pressure_nn)./delta_xcb;
anv.visc(1:Nx) = 0.0; %viscosity*(2.0d+0 * neutral_velocity_cb(0) + neutral_velocity(2)-3.0d+0*neutral_velocity(1))  
anv.visc(2:Nx-1) = in.numerics.viscosity*(o.neutral_velocity(Ntime,3:Nx) + o.neutral_velocity(Ntime,1:Nx-2) - 2*o.neutral_velocity(Ntime,2:Nx-1));
%ydot(4*Nx+ix) + viscosity*(neutral_velocity(ix+1) + neutral_velocity(ix-1)-2.0d+0*neutral_velocity(ix))
anv.src = o.Source_vn(Ntime,:);
anv.dt = anv.dgamdx +anv.dpdx + anv.src+ anv.visc;
error.bal.neutral_mom = anv.dt./pressure_nn_cc;
error.src.neutral_mom = anv.src*nan;

% 6. molecules
mol = struct;
mol.dgamdx = -diff(o.Gamma_molecule(Ntime,:))./delta_xcb;
mol.src = o.Source_molecule(Ntime,:);
mol.dt = mol.dgamdx + mol.src;
error.bal.molecule = mol.dt ./o.molecule(Ntime,:);
error.src.molecule = mol.src*nan; % this should be completed

% 7. core particle balance
cor = struct;
cor.gam2sol = -o.Gamma_core2sol(Ntime,1);
cor.ext2cora = o.sum_extern2core_flux(Ntime,1);
cor.ext2corm = o.sum_extern2core_mol(Ntime,1);
cor.sol2cora = o.sol2core_flux(Ntime,1);
cor.sol2corm = o.sol2core_mol(Ntime,1);
cor.fuelling = in.dynamic.dyn_core_fuelling(Ntime,1);
cor.dt = +cor.gam2sol + cor.ext2cora + cor.ext2corm + cor.sol2cora + cor.sol2corm+cor.fuelling;
cor.error = cor.dt./o.core_density(Ntime,1);

% 8. fluxes traded between domains 
trade = struct;
trade.coresol.c2s_flux = sum(o.core2sol_flux(Ntime,:));
trade.coresol.c2s_mol = sum(o.core2sol_mol(Ntime,:));
trade.coresol.s2c_flux = o.sol2core_flux(Ntime,1);
trade.coresol.s2c_mol =  o.sol2core_mol(Ntime,1);
trade.coresol.sum = trade.coresol.c2s_flux + trade.coresol.c2s_mol + trade.coresol.s2c_flux + trade.coresol.s2c_mol;
trade.coresol.norm = max(abs([trade.coresol.c2s_flux, trade.coresol.c2s_mol , trade.coresol.s2c_flux , trade.coresol.s2c_mol]));
trade.error.coresol = max(abs(trade.coresol.sum./trade.coresol.norm),0);

trade.coreext.e2c_flux = -sum(o.extern2core_flux(Ntime,:));
trade.coreext.e2c_mol = -sum(o.extern2core_mol(Ntime,:));
trade.coreext.se2c_flux = sum(o.sum_extern2core_flux(Ntime,1));
trade.coreext.se2c_mol = sum(o.sum_extern2core_mol(Ntime,1));
trade.coreext.sum = trade.coreext.e2c_flux+trade.coreext.e2c_mol + trade.coreext.se2c_flux + trade.coreext.se2c_mol;
trade.coreext.norm =  max(abs([trade.coreext.e2c_flux,trade.coreext.e2c_mol ,trade.coreext.se2c_flux , trade.coreext.se2c_mol]));
trade.error.coreext = max(abs(trade.coreext.sum./trade.coreext.norm),0);

trade.extsol.e2s_flux = sum(o.extern2sol_flux(Ntime,:));
trade.extsol.e2s_mol =  sum(o.extern2sol_mol(Ntime,:));
trade.extsol.s2e_flux = sum(o.sol2extern_flux(Ntime,:));
trade.extsol.s2e_mol = sum(o.sol2extern_mol(Ntime,:));
trade.extsol.sum = trade.extsol.s2e_flux + trade.extsol.s2e_mol + trade.extsol.e2s_flux + trade.extsol.e2s_mol;
trade.extsol.norm = max(abs([trade.extsol.s2e_flux , trade.extsol.s2e_mol , trade.extsol.e2s_flux , trade.extsol.e2s_mol]));
trade.error.extsol = max(abs(trade.extsol.sum./trade.extsol.norm),0);

% reserrvoir balances
% flows between neutral chambers
        exaflow = o.extern_neutral_flux(Ntime,:);
         %!extern_fluxes =          (/(extern_neutral(5)-extern_neutral(1))*extern_ex(1), &   ! from 5->1
         %!                         (extern_neutral(2)-extern_neutral(3))*extern_ex(2), &	    ! from 2->3
         %!                         (extern_neutral(3)-extern_neutral(4))*extern_ex(3) /)     ! from 3->4
         	%(/ extern_neutral_flux(1) , & ! positive direction for neutral flows outside the plasma is anti-clockwise
			%  -extern_neutral_flux(2) , &
			%   extern_neutral_flux(2) - extern_neutral_flux(3) , &
			 %  extern_neutral_flux(3) , &
			 % !-extern_neutral_flux(1) /)
        ext_flows_ato(1:5) = [exaflow(1), ...
                                    -exaflow(2),  ...
                                     exaflow(2)-exaflow(3), ...
                                     exaflow(3), ...
                                    -exaflow(1)];
        exmflow = o.extern_molecule_flux(Ntime,:);
        ext_flows_mol(1:5) = [exmflow(1), ...
                                    -exmflow(2),  ...
                                     exmflow(2)-exmflow(3), ...
                                     exmflow(3), ...
                                    -exmflow(1)];

reservoir.ato.pump = -o.extern_neutral_density(Ntime,:).*in.physics.pump_rate_n.*in.physics.extern_neutral_volumes;
reservoir.ato.ext2ext_flux = ext_flows_ato(1:5) ;
reservoir.ato.sol2ext_flux = o.sol2extern_flux(Ntime,:);
reservoir.ato.tar2ext_flux = o.tar2extern_flux(Ntime,:);
reservoir.ato.ext2cor_flux = -o.extern2core_flux(Ntime,:);
reservoir.ato.puff = in.dynamic.dyn_neutral_puff(Ntime,:) ;
reservoir.ato.dt = reservoir.ato.sol2ext_flux ...
                 + reservoir.ato.tar2ext_flux ...
                 + reservoir.ato.ext2cor_flux ... 
                 + in.dynamic.dyn_neutral_puff(Ntime,:) ...
                 + reservoir.ato.pump ...
                 + reservoir.ato.ext2ext_flux;
reservoir.ato.error = reservoir.ato.dt./o.extern_neutral_density(Ntime,:);

reservoir.mol.pump = -o.extern_molecule_density(Ntime,:).*in.physics.pump_rate_m.*in.physics.extern_neutral_volumes;
reservoir.mol.ext2ext_mol = ext_flows_mol;
reservoir.mol.sol2ext_mol = o.sol2extern_mol(Ntime,:);
reservoir.mol.tar2ext_mol = o.tar2extern_mol(Ntime,:);
reservoir.mol.ext2cor_mol = -o.extern2core_mol(Ntime,:);
reservoir.mol.puff = in.dynamic.dyn_molecule_puff(Ntime,:) ;
reservoir.mol.dt = reservoir.mol.sol2ext_mol...
                 + reservoir.mol.tar2ext_mol ...
                 + reservoir.mol.ext2cor_mol ...
                 + reservoir.mol.puff ...
                 + reservoir.mol.pump...
                 + reservoir.mol.ext2ext_mol;
reservoir.mol.error = reservoir.mol.dt./o.extern_molecule_density(Ntime,:);

% sol reservoir balance
reservoir.sol.ext2sol_flux = sum(o.extern2sol_flux(Ntime,:));
reservoir.sol.ext2sol_mol = 2*sum(o.extern2sol_mol(Ntime,:));
reservoir.sol.cor2sol_flux = sum(o.core2sol_flux(Ntime,:));
reservoir.sol.cor2sol_mol = 2*sum(o.core2sol_mol(Ntime,:));
reservoir.sol.tar2ext_flux = -sum(o.tar2extern_flux(Ntime,:));
reservoir.sol.tar2ext_mol = -2*sum(o.tar2extern_mol(Ntime,:));
reservoir.sol.farsol_ion = o.sum_sol2extern_ion_mol(Ntime) + o.sum_sol2extern_ion_flux(Ntime);
reservoir.sol.gamma_core = in.dynamic.dyn_gam_core(Ntime);
reservoir.sol.dt = reservoir.sol.ext2sol_mol ...
                + reservoir.sol.ext2sol_flux...
                + reservoir.sol.cor2sol_flux ...
                + reservoir.sol.cor2sol_mol ...
                + reservoir.sol.tar2ext_flux...
                + reservoir.sol.tar2ext_mol...
                + reservoir.sol.farsol_ion...
                + reservoir.sol.gamma_core;


% balances
bal = struct;
bal.den = den;
bal.mom = mom;
bal.ene = ene;
bal.ato = ato;
bal.anv = anv;
bal.mol = mol;


 % return processed terms    
    proc.pressure = pressure_ss;
    proc.P_kin = P_kin;
    proc.P_tot = P_tot;  
    proc.pressure_nn = pressure_nn;
    % post process the AMJUEL rates
    % --    todo --

    Nx = length(x);

fprintf('Note: the following error is between DIV1D and what we would expect \n      from a parallel calculation in Matlab (which is slightly different).\n\n')
fprintf('Average normalized errors \n')
fprintf('density  :  %f \n',sum(error.bal.density)/Nx);
fprintf('momentum :  %f \n',sum(error.bal.momentum)/Nx);
fprintf('energy   :  %f \n',sum(error.bal.energy)/Nx);
fprintf('neutral  :  %f \n',sum(error.bal.neutral)/Nx);
fprintf('exchange : c2s %e, c2e %e, e2s %e \n', trade.error.coresol,trade.error.coreext,trade.error.extsol)
% if(in.numerics.evolve_core == 1) % perhaps turn off iff not applicable
fprintf('core     : %e \n', cor.error);
% end
fprintf('ato res  : %f %f %f %f %f \n',reservoir.ato.error);
fprintf('mol res  : %f %f %f %f %f \n \n',reservoir.mol.error);

%% check boundary conditions
% the boundary conditions in div1d are set for the fluxes at the target:
%   density flux: gamma_n = density*velocity_target 
%  momentum flux: gamma_mom = density*mass*velocity**2
%    energy flux: q_parallel = gamma*velocity*density*e_charge*T
% In the following we compare:
%   1. the output of DIV1D for these fluxes at the last cell boundary (i.e. target)
%   2. the fluxes that should be enforced in DIV1D based on n,T at cell centers.
%   3. fluxes as calculated from extrapolating n,T,v to the last cell boundary
   
% >>> 2. fluxes as calculated in DIV1D with equidistant discretization
csound_target_c = sqrt(2.0*e_charge*(1.5*o.temperature(Ntime,end)-0.5*o.temperature(Ntime,end-1))/in.physics.mass);
velocity_c      = max(o.velocity(Ntime,end),csound_target_c);
density_c       = min(o.density(Ntime,end),(1.5*o.density(Ntime,end)-0.5*o.density(Ntime,end-1)));
temperature_c   = max(1.5*o.temperature(Ntime,end)-0.5*o.temperature(Ntime,end-1),0.1);
gamma_target_c  = density_c*velocity_c;
gamma_mom_c     = density_c*in.physics.mass*velocity_c^2;
q_target_c        = in.physics.gamma*csound_target_c*density_c*e_charge*temperature_c;

% >>> 3. fluxes as calculated in DIV1D with non-equidistant discretization (from version 5.0 onward)
csound_target_d = sqrt(2.0*e_charge*(temperature_cb(end)/in.physics.mass));
velocity_d      = max(velocity_cb(end),csound_target_d);
gamma_target_d  = density_cb(end)*velocity_d;
gamma_mom_d     = density_cb(end)*in.physics.mass*velocity_d^2;
q_target_d      = in.physics.gamma*csound_target_d(end)*density_cb(end)*e_charge*temperature_cb(end);

fprintf('Note: the following boundary conditions checks are on code implementation \n      and values calculated from extrapolating with Matlab \n \n')

fprintf('BC checks: DIV1D-output  check-eqidist  check-non-equidist    \n')
fprintf('gamma_n :  %e  %e  %e \n',o.Gamma_n(Ntime,end), gamma_target_c, gamma_target_d)
fprintf('gamma_mom: %e  %e  %e \n',o.Gamma_mom(Ntime,end),gamma_mom_c,gamma_mom_d)
fprintf('q_parallel:%e  %e  %e \n \n',o.q_parallel(Ntime,end),q_target_c,q_target_d)
fprintf('Bohm condition \n')
fprintf('DIV1D %e, equidist %e, non-equidist %e \n \n',velocity_cb(end),csound_target_c,csound_target_d)


%% do 2PM calculations for DIV1D to interpret outputs and check for code errors
% determine geometry
% !!! this should be updated!!! 
if in.physics.l_core_sol > 0.0 
    if in.physics.x_core_sol > 0.0 %  two targets
    geometry = 'double';
    xind = in.grid.i_xpoint(2);
    q_u_t = o.q_parallel(:,1)*0 +1
    % curious on how to deal with the signs here
    % forget about the inner target for now.
    else % one target with core-sol
    geometry = 'singlecore';
    xind = in.grid.i_xpoint(2);
    q_u_t = o.q_parallel(:,1)*0 +1; % how to deal with this... 
    end
else % one target no core
geometry = 'single';
xind = 1;
q_u_t = o.q_parallel(:,xind);
end
q_X_t = q_u_t; %q_parallel(:,xind);

switch geometry
    case 'single'
    q_X = q_X_t(Ntime);
    for Ntime = 1: length(o.time)
    D2PM.n_X(Ntime) = o.density(Ntime,1);
    D2PM.n_L(Ntime) = mean(o.density(Ntime,Nx-30:Nx)); 
    D2PM.T_L(Ntime) = 1.5*o.temperature(Ntime,Nx)-0.5*o.temperature(Ntime,Nx-1);
    D2PM.T_X(Ntime) = o.temperature(Ntime,1);
    D2PM.q_X(Ntime) = abs(q_X);
    D2PM.L   = L;
    % convective fraction
    f_cnv = (5*D2PM.n_X(Ntime)*e_charge*D2PM.T_X(Ntime)*o.velocity(Ntime,1))/D2PM.q_X(Ntime);
    D2PM.f_cnv(Ntime) = max(min(0.9999,f_cnv),0);
    
    % power and momentum loss fractions with volume integral
    f_pwr = abs(-sum(delta_xcb.*Source_Qtot)/q_X);
    f_mom = abs(-sum(delta_xcb.*(1./B_field*m.*o.density(Ntime,:).*o.velocity(Ntime,:).^2.*[diff(B_field) 0 ]./delta_xcb ...
                                    + o.Source_v(Ntime,:)...
                                    + [0 in.numerics.viscosity*...
      (o.velocity(Ntime,3:end)-2*o.velocity(Ntime,2:end-1)+o.velocity(Ntime,1:end-2)) 0] ))...
                   /total_pressure(Ntime,1));      
    D2PMwivi.f_mom(Ntime)=max(min(0.9999,f_mom),0);
    D2PMwivi.f_pwr(Ntime)=max(min(0.9999,f_pwr),0);
    
    str = '';% predictions by KR and 2PM Westerhof
    tmp = test_2PM_KR(D2PMwivi.f_pwr(Ntime),D2PMwivi.f_mom(Ntime),D2PM.n_X(Ntime),D2PM.n_L(Ntime), ...
                      D2PM.T_X(Ntime),D2PM.T_L(Ntime),D2PM.q_X(Ntime),D2PM.L, ...
                      'print',str,'gamma',P.gamma,'f_exp',flux_expansion);
    D2PMwivi.T_X(Ntime,1:3) = tmp.T_u; 
    D2PMwivi.T_L(Ntime,1:3) = tmp.T_t; 
    D2PMwivi.n_L(Ntime,1:3) = tmp.n_t;   
    
    % power and momentum loss fraction without volume integral
    f_mom(Ntime) = abs(total_pressure(Ntime,1)-mean(total_pressure(Ntime,Nx-5:Nx)))./total_pressure(Ntime,1);
    D2PMwovi.f_mom(Ntime)=max(min(0.9999,f_mom(Ntime)),0);
    % fpwr = 1 - T_t * Gam_t / (qu*e*gamma)  * Rt/Ru  = 1 - T_t * Gam_t / (qu*e*gamma) * epsilon_fluxexpansion
    
    %q_plasma_cooling  =  gammasolps*k*T_t Gamma_t
    q_plasma_cool       = P.gamma*e_charge* D2PM.T_L(Ntime)*mean(o.Gamma_n(Ntime,Nx-10:Nx)); D2PM.q_L(Ntime) = q_plasma_cool;
    % fpwr = 1 - q_plasma_cool/qu * rt/ru   % ( eq 64 in Stangeby 2018 )
    f_pwr(Ntime) = 1 - q_plasma_cool/abs(D2PM.q_X(Ntime)) * FE2PM; % *(1/S_flux_expansion);
    
    % f_pwr(Ntime) = 1 - mean(div1doutput.Gamma_n(Ntime,Nx-10:Nx)) * D2PM.T_L(Ntime)*e_charge*div1dinput.gamma/(abs(D2PM.q_X(Ntime))/flux_expansion); % Ru/Rt    % (eq 43 Stangeby 2018) 
    D2PMwovi.f_pwr(Ntime) = max(min(f_pwr(Ntime),0.9999),0);  
    
    str = '';% predictions by KR and 2PM westerhof
    tmp = test_2PM_KR(D2PMwovi.f_pwr(Ntime),D2PMwovi.f_mom(Ntime),D2PM.n_X(Ntime),D2PM.n_L(Ntime), ...
                      D2PM.T_X(Ntime),D2PM.T_L(Ntime),D2PM.q_X(Ntime),D2PM.L, ...
                      'print',str,'gamma',P.gamma,'f_exp',FE2PM);
    D2PMwovi.T_X(Ntime,1:3) = tmp.T_u; 
    D2PMwovi.T_L(Ntime,1:3) = tmp.T_t; 
    D2PMwovi.n_L(Ntime,1:3) = tmp.n_t;   
    end
    
    Ntime = P.Ntime; % print the predictions at the requested timesteps
    if P.print ==1; str='DIV1D with volume integral';
    [~] = test_2PM_KR(D2PMwivi.f_pwr(Ntime),D2PMwivi.f_mom(Ntime),D2PM.n_X(Ntime),D2PM.n_L(Ntime), ...
                      D2PM.T_X(Ntime),D2PM.T_L(Ntime),D2PM.q_X(Ntime),D2PM.L, ...
                      'print',str,'gamma',P.gamma,'f_cnv',D2PM.f_cnv(Ntime),'f_exp',FE2PM);              
    end
    if P.print ==1; str='DIV1D without volume integral';
    [~] = test_2PM_KR(D2PMwovi.f_pwr(Ntime),D2PMwovi.f_mom(Ntime),D2PM.n_X(Ntime),D2PM.n_L(Ntime), ...
                      D2PM.T_X(Ntime),D2PM.T_L(Ntime),D2PM.q_X(Ntime),D2PM.L, ...
                      'print',str,'gamma',P.gamma,'f_cnv',D2PM.f_cnv(Ntime),'f_exp',FE2PM);              
    end

    case 'singlecore'
    xind = in.grid.i_xpoint(2);
    q_X = q_X_t(Ntime); 
    for Ntime = 1: length(o.time)
    D2PM.n_X(Ntime) = o.density(Ntime,xind);
    D2PM.n_L(Ntime) = mean(o.density(Ntime,Nx-30:Nx)); 
    D2PM.T_L(Ntime) = 1.5*o.temperature(Ntime,Nx)-0.5*o.temperature(Ntime,Nx-1);
    D2PM.T_X(Ntime) = o.temperature(Ntime,xind);
    D2PM.q_X(Ntime) = abs(q_X);
    D2PM.L   = L-L_core_SOL;
    % convective fraction
    f_cnv = (5*D2PM.n_X(Ntime)*e_charge*D2PM.T_X(Ntime)*o.velocity(Ntime,xind))/D2PM.q_X(Ntime);
    D2PM.f_cnv(Ntime) = max(min(0.9999,f_cnv),0);
    
    % power and momentum loss fractions with volume integral
    f_pwr = abs(-sum(delta_xcb(xind:end).*Source_Qtot(xind:end))/q_X);
    f_mom = abs(-sum(delta_xcb(xind:end).*(1./B_field(xind:end)*m.*o.density(Ntime,xind:end).*o.velocity(Ntime,xind:end).^2.*[diff(B_field_cb(xind:end))]./delta_xcb(xind:end) ...
                                    + o.Source_v(Ntime,xind:end)...
                                    + [in.numerics.viscosity*...
      (o.velocity(Ntime,xind+1:end)-2*o.velocity(Ntime,xind:end-1)+o.velocity(Ntime,xind-1:end-2)) 0] ))...
                   /total_pressure(Ntime,xind));      
    D2PMwivi.f_mom(Ntime)=max(min(0.9999,f_mom),0);
    D2PMwivi.f_pwr(Ntime)=max(min(0.9999,f_pwr),0);
    
    str = '';% predictions by KR and 2PM Westerhof
    tmp = test_2PM_KR(D2PMwivi.f_pwr(Ntime),D2PMwivi.f_mom(Ntime),D2PM.n_X(Ntime),D2PM.n_L(Ntime), ...
                      D2PM.T_X(Ntime),D2PM.T_L(Ntime),D2PM.q_X(Ntime),D2PM.L, ...
                      'print',str,'gamma',P.gamma,'f_exp',flux_expansion);
    D2PMwivi.T_X(Ntime,1:3) = tmp.T_u; 
    D2PMwivi.T_L(Ntime,1:3) = tmp.T_t; 
    D2PMwivi.n_L(Ntime,1:3) = tmp.n_t;   
    
    % power and momentum loss fraction without volume integral
    f_mom(Ntime) = abs(total_pressure(Ntime,xind)-mean(total_pressure(Ntime,Nx-5:Nx)))./total_pressure(Ntime,xind);
    D2PMwovi.f_mom(Ntime)=max(min(0.9999,f_mom(Ntime)),0);
    % fpwr = 1 - T_t * Gam_t / (qu*e*gamma)  * Rt/Ru  = 1 - T_t * Gam_t / (qu*e*gamma) * epsilon_fluxexpansion
    
    %q_plasma_cooling  =  gammasolps*k*T_t Gamma_t
    q_plasma_cool       = P.gamma*e_charge* D2PM.T_L(Ntime)*mean(o.Gamma_n(Ntime,Nx-10:Nx)); D2PM.q_L(Ntime) = q_plasma_cool;
    % fpwr = 1 - q_plasma_cool/qu * rt/ru   % ( eq 64 in Stangeby 2018 )
    f_pwr(Ntime) = 1 - q_plasma_cool/abs(D2PM.q_X(Ntime)) * FE2PM; % *(1/S_flux_expansion);
    
    % f_pwr(Ntime) = 1 - mean(div1doutput.Gamma_n(Ntime,Nx-10:Nx)) * D2PM.T_L(Ntime)*e_charge*div1dinput.gamma/(abs(D2PM.q_X(Ntime))/flux_expansion); % Ru/Rt    % (eq 43 Stangeby 2018) 
    D2PMwovi.f_pwr(Ntime) = max(min(f_pwr(Ntime),0.9999),0);  
    
    str = '';% predictions by KR and 2PM westerhof
    tmp = test_2PM_KR(D2PMwovi.f_pwr(Ntime),D2PMwovi.f_mom(Ntime),D2PM.n_X(Ntime),D2PM.n_L(Ntime), ...
                      D2PM.T_X(Ntime),D2PM.T_L(Ntime),D2PM.q_X(Ntime),D2PM.L, ...
                      'print',str,'gamma',P.gamma,'f_exp',FE2PM);
    D2PMwovi.T_X(Ntime,1:3) = tmp.T_u; 
    D2PMwovi.T_L(Ntime,1:3) = tmp.T_t; 
    D2PMwovi.n_L(Ntime,1:3) = tmp.n_t;   
    end
    
    Ntime = P.Ntime; % print the predictions at the requested timesteps
    if P.print ==1; str='DIV1D with volume integral';
    [~] = test_2PM_KR(D2PMwivi.f_pwr(Ntime),D2PMwivi.f_mom(Ntime),D2PM.n_X(Ntime),D2PM.n_L(Ntime), ...
                      D2PM.T_X(Ntime),D2PM.T_L(Ntime),D2PM.q_X(Ntime),D2PM.L, ...
                      'print',str,'gamma',P.gamma,'f_cnv',D2PM.f_cnv(Ntime),'f_exp',FE2PM);              
    end
    if P.print ==1; str='DIV1D without volume integral';
    [~] = test_2PM_KR(D2PMwovi.f_pwr(Ntime),D2PMwovi.f_mom(Ntime),D2PM.n_X(Ntime),D2PM.n_L(Ntime), ...
                      D2PM.T_X(Ntime),D2PM.T_L(Ntime),D2PM.q_X(Ntime),D2PM.L, ...
                      'print',str,'gamma',P.gamma,'f_cnv',D2PM.f_cnv(Ntime),'f_exp',FE2PM);              
    end
    case 'double'
    % not implemented yet

end

%% Plotting
% if P.plotoutput
%  plot_div1d_v600(o,in)
% end
if P.plotbal 
    plot_div1d_bal_v600(o.X,bal);
end
if P.ploterror
 plot_div1d_error_v600(o.X,error);
end
if P.plotprofiles
  plot_div1d_profiles_v600(o, in,'avg',1,'plot',1,'figtight',1,'Ntime',P.Ntime,'print',1,...  );
     'save',0,'allpossiblefigures',0,'hold',1,'plot2PM',0,'xlimits',P.xlimits,'version',P.version,'FE',in.physics.flux_expansion,'cratio',P.cratio)
end
% if P.plotreservoirs
%     plot_div1d_reservoirs_v600(o,in);
% end

%% output
out = struct;
out.proc = proc;
out.bal = bal; % balances
out.cor = cor; % core balance
out.trade = trade; % the trading fluxes should all be zero...
out.reservoir = reservoir; % the checks on reservoirs at Ntime
try
out.D2PM = D2PM;
out.D2PMwivi = D2PMwivi; % with volumetric integration
out.D2PMwovi = D2PMwovi; % without volume integration
catch
 disp('failed to do 2P analysis')
end
out.Derror = error;

end

