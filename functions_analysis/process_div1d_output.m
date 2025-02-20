function [out] = process_div1d_output(output,input,varargin)
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
D.gamma = input.physics.gamma;
D.cratio = 0.03;
% -- numerics ---
%D.gamma = input.numerics.6;
D.dxmin = 0.05; 
D.xlimits = [0  7.8];
D.Ntime = length(output.time);
D.solps = 0;
D.input = 0;
D.version = 'v2.0.0';

% -- plotting ---
D.plotoutput = 0;
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
Nx          = length(output.X);
x           = output.X;
xcb         = output.Xcb;
L           = input.physics.l;
L_core_SOL  = input.physics.l_core_sol;
x_core_SOL  = input.physics.x_core_sol;

B_field = input.grid.b_field;
B_field_cb = input.grid.b_field_cb;
B_tot = input.grid.b_field.*input.grid.b_trans;
B_tot_cb = input.grid.b_field_cb.*input.grid.b_trans_cb;

P.tn            = input.physics.neutral_residence_time;
flux_expansion  = input.physics.flux_expansion;
if P.FE ==1
    FE2PM = flux_expansion;
    P.FE = input.physics.flux_expansion;
else
    FE2PM = P.FE; % flux expansion factor used in Two-Point Formatting
end

% calculate spacing of the grid
delta_x = diff(x); 
delta_xcb = diff(xcb);

% extract important profiles
density         = output.density;
density_cb  = interp1(x,density(Ntime,:),xcb,"linear", "extrap" );
temperature     = output.temperature;
temperature_cb  = interp1(x,temperature(Ntime,:),xcb,"linear", "extrap" );
q_parallel      = output.q_parallel;
neutral_density = output.neutral_density;
Source_v        = output.Source_v;
Gamma_mom       = output.Gamma_mom;
Gamma_n         = output.Gamma_n;
neutral_flux    = output.neutral_flux;
velocity        = output.velocity;
velocity_cb     = interp1(x,velocity(Ntime,:),xcb,"linear", "extrap" ); 
Source_Q        = output.Source_Q;
Source_n        = output.Source_n;
Source_neutral  = output.Source_neutral;

% momentum
total_pressure = 2.0*density.*temperature*e_charge + (Gamma_mom(:,2:end)+Gamma_mom(:,1:end-1))/2;
pressure_t     = 2.0*density.*temperature*e_charge;
pressure_ss    = pressure_t(Ntime,:);
P_kin = 2.0*density(Ntime,:).*temperature(Ntime,:)*e_charge;
P_ram = (Gamma_mom(Ntime,2:end)+Gamma_mom(Ntime,1:end-1))/2; %.*B_field
P_tot = P_kin + P_ram;

% energy
Source_Qtot = Source_Q(Ntime,:);
Source_Qtot(2:Nx) = Source_Qtot(2:Nx) + velocity(Ntime,2:Nx).*(pressure_ss(2:Nx)-pressure_ss(1:Nx-1))./delta_x(1:Nx-1);
Source_Qtot = Source_Qtot./B_tot; 


%% post process rates
for i = length(output.X):-1:1
    nenn = density(Ntime,i).*neutral_density(Ntime,i);
    rate_exc(i)     = nenn*AMJUEL_H102_R215_e_loss_excitation(density(Ntime,i),temperature(Ntime,i)).*e_charge; 
    rate_rec(i)     = density(Ntime,i)*density(Ntime,i)*AMJUEL_H46_R218_recombination(density(Ntime,i),temperature(Ntime,i));
    rate_ion_vel(i) = nenn*AMJUEL_H43_R215_ionization(density(Ntime,i),temperature(Ntime,i)).*0.5*input.physics.mass*velocity(Ntime,i)^2;
    rate_ion_nen(i) = nenn*AMJUEL_H43_R215_ionization(density(Ntime,i),temperature(Ntime,i)).*e_charge.*abs(input.physics.neutral_energy).*temperature(Ntime,i);
    rate_ion_den(i) = nenn*AMJUEL_H43_R215_ionization(density(Ntime,i),temperature(Ntime,i));
    rate_cx(i)      = nenn*AMJUEL_H219_R318_charge_exchange(temperature(Ntime,i),input.physics.mass); % charge exchange requires mass to scale the rate.
    rate_rec_ene(i)    = density(Ntime,i)*density(Ntime,i)*AMJUEL_104_218_e_loss_recombination(density(Ntime,i), temperature(Ntime,i));
    rate_cooling_post_c = density(Ntime,i)*density(Ntime,i)*input.physics.impurity_concentration(1)*POST_cooling_carbon(temperature(Ntime,i));
end
    proc.rate_cx        = rate_cx;
    proc.rate_exc       = rate_exc;
    proc.rate_rec       = rate_rec;
    proc.rate_ion_vel   = rate_ion_vel;
    proc.rate_ion_nen   = rate_ion_nen;
    proc.rate_ion_den   = rate_ion_den;
    proc.rate_rec_ene   = rate_rec_ene;
    proc.rate_cooling_post_c = rate_cooling_post_c;

%% post process sources
% TODO: check the source terms with these rates

%% Perform checks on stationarity and accuracy (DIV1D)
% 1. density equation
    error = struct;
    den_rhs = -B_tot.*diff(Gamma_n(Ntime,:)./B_tot_cb)./delta_xcb; % -B d/dx (Gamma/B )
    error.dndt = (den_rhs + Source_n(Ntime,:));
    error.density = error.dndt./density(Ntime,:);
    error.source_n = Source_n(Ntime,:) - rate_ion_den + rate_rec;
% 2. momentum equation
    mom_bdgamdx = -B_tot.*diff(Gamma_mom(Ntime,:)./B_tot_cb)./delta_xcb;
    pressure_cb = interp1(x,pressure_ss,xcb,'linear','extrap');
    mom_dpdx = -diff(pressure_cb)./delta_xcb;
    mom_visc = zeros(1,Nx);
    mom_visc(2:Nx-1) = input.numerics.viscosity*(velocity(Ntime,3:end) + velocity(Ntime,1:end-2) - 2*velocity(Ntime,2:end-1));
    error.dmomdt = mom_dpdx + mom_bdgamdx + Source_v(Ntime,:);
    error.momentum = error.dmomdt./pressure_ss; 
    error.source_v = Source_v(Ntime,:) + input.physics.mass*velocity(Ntime,1:end).*(rate_cx - rate_rec) ;
% % 3. energy equation
    ene_bdqdx = -B_tot.*diff(q_parallel(Ntime,:)./B_tot_cb)./delta_xcb;
    velocity_cc = interp1(x,velocity(Ntime,:),x);
    ene_vdpdx = velocity_cc .* diff(pressure_cb)./delta_xcb;
    error.dqdt = ene_bdqdx + ene_vdpdx  + Source_Q(Ntime,:);
    error.energy = error.dqdt./(3.0.*density(Ntime,:).*e_charge.*temperature(Ntime,:));
    error.source_Q = Source_Q(Ntime,:)- ( -1.5*e_charge.*temperature(Ntime,:).*2.*rate_rec  ...
                     -(1.5+input.physics.neutral_energy)*e_charge.*temperature(Ntime,:).*rate_cx ...
                     +rate_ion_nen ...
                     +rate_ion_vel ...
                     -rate_exc ...
                     +(rate_rec*13.6 - rate_rec_ene).*e_charge ...
                     - rate_cooling_post_c+ ...
                     + input.grid.core_source_profile*input.physics.q_core./input.grid.volumes);
      %error.source_Q = Source_Q(Ntime,:) + 1.5*e_charge*temperature(Ntime,:).*(rate_cx +2*rate_rec) + ...
      %               -rate_exc +rate_rec*12.6*e_charge - rate_rec_ene - rate_cooling_post_c;
% 4. neutral equation
    neu_diffdnndx = -diff(neutral_flux(Ntime,:))./delta_xcb;
    error.dnndt = neu_diffdnndx + Source_neutral(Ntime,:);
    error.neutral = error.dnndt./neutral_density(Ntime,:);
    error.Source_neutral = Source_neutral(Ntime,:) + Source_n(Ntime,:);

 % return processed terms    
    proc.pressure = pressure_ss;
    proc.P_kin = P_kin;
    proc.P_tot = P_tot;  
    proc.ene_bdqdx = ene_bdqdx;
    proc.ene_vdpdx = ene_vdpdx;
    proc.mom_dpdx = mom_dpdx;
    proc.mom_visc = mom_visc;
    proc.mom_bdgamdx = mom_bdgamdx;
    proc.den_rhs= den_rhs;
    proc.neu_diffdnndx = neu_diffdnndx;
    proc.dqdt = error.dqdt;
    proc.dmomdt = error.dmomdt;
    proc.dndt = error.dndt;
    proc.dnndt = error.dnndt;
    % post process the AMJUEL rates
    % --    todo --

    Nx = length(x);

fprintf('Note: the following error is between DIV1D and what we would expect \n      from a parallel calculation in Matlab (which is slightly different).\n\n')
fprintf('Average normalized errors \n')
fprintf('density  :  %f \n',sum(error.density)/Nx);
fprintf('momentum :  %f \n',sum(error.momentum)/Nx);
fprintf('energy   :  %f \n',sum(error.energy)/Nx);
fprintf('neutral  :  %f \n \n',sum(error.neutral)/Nx);

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
csound_target_c = sqrt(2.0*e_charge*(1.5*temperature(Ntime,end)-0.5*temperature(Ntime,end-1))/input.physics.mass);
velocity_c      = max(velocity(Ntime,end),csound_target_c);
density_c       = min(density(Ntime,end),(1.5*density(Ntime,end)-0.5*density(Ntime,end-1)));
temperature_c   = max(1.5*temperature(Ntime,end)-0.5*temperature(Ntime,end-1),0.1);
gamma_target_c  = density_c*velocity_c;
gamma_mom_c     = density_c*input.physics.mass*velocity_c^2;
q_target_c        = input.physics.gamma*csound_target_c*density_c*e_charge*temperature_c;

% >>> 3. fluxes as calculated in DIV1D with non-equidistant discretization (from version 5.0 onward)
csound_target_d = sqrt(2.0*e_charge*(temperature_cb(end)/input.physics.mass));
velocity_d      = max(velocity_cb(end),csound_target_d);
gamma_target_d  = density_cb(end)*velocity_d;
gamma_mom_d     = density_cb(end)*input.physics.mass*velocity_d^2;
q_target_d      = input.physics.gamma*csound_target_d(end)*density_cb(end)*e_charge*temperature_cb(end);

fprintf('Note: the following boundary conditions checks are on code implementation \n      and values calculated from extrapolating with Matlab \n \n')

fprintf('BC checks: DIV1D-output  check-eqidist  check-non-equidist    \n')
fprintf('gamma_n :  %e  %e  %e \n',Gamma_n(Ntime,end), gamma_target_c, gamma_target_d)
fprintf('gamma_mom: %e  %e  %e \n',Gamma_mom(Ntime,end),gamma_mom_c,gamma_mom_d)
fprintf('q_parallel:%e  %e  %e \n \n',q_parallel(Ntime,end),q_target_c,q_target_d)
fprintf('Bohm condition \n')
fprintf('DIV1D %e, equidist %e, non-equidist %e \n \n',velocity_cb(end),csound_target_c,csound_target_d)


%% do 2PM calculations for DIV1D to interpret outputs and check for code errors
% determine geometry
% !!! this should be updated!!! 
if input.physics.l_core_sol > 0.0 
    if input.physics.x_core_sol > 0.0 %  two targets
    geometry = 'double';
    xind = input.grid.i_xpoint;
    % curious on how to deal with the signs here
    else % one target with core-sol
    geometry = 'singlecore';
    xind = input.grid.i_xpoint(2);
    q_u_t = q_parallel(:,1)*0 +1; % how to deal with this... 
    end
else % one target no core
geometry = 'single';
xind = 1;
q_u_t = q_parallel(:,xind);
end
q_X_t = q_u_t; %q_parallel(:,xind);

switch geometry
    case 'single'
    q_X = q_X_t(Ntime);
    for Ntime = 1: length(output.time)
    D2PM.n_X(Ntime) = output.density(Ntime,1);
    D2PM.n_L(Ntime) = mean(output.density(Ntime,Nx-30:Nx)); 
    D2PM.T_L(Ntime) = 1.5*output.temperature(Ntime,Nx)-0.5*output.temperature(Ntime,Nx-1);
    D2PM.T_X(Ntime) = output.temperature(Ntime,1);
    D2PM.q_X(Ntime) = abs(q_X);
    D2PM.L   = L;
    % convective fraction
    f_cnv = (5*D2PM.n_X(Ntime)*e_charge*D2PM.T_X(Ntime)*output.velocity(Ntime,1))/D2PM.q_X(Ntime);
    D2PM.f_cnv(Ntime) = max(min(0.9999,f_cnv),0);
    
    % power and momentum loss fractions with volume integral
    f_pwr = abs(-sum(delta_xcb.*Source_Qtot)/q_X);
    f_mom = abs(-sum(delta_xcb.*(1./B_field*m.*density(Ntime,:).*velocity(Ntime,:).^2.*[diff(B_field) 0 ]./delta_xcb ...
                                    + Source_v(Ntime,:)...
                                    + [0 input.numerics.viscosity*...
      (velocity(Ntime,3:end)-2*velocity(Ntime,2:end-1)+velocity(Ntime,1:end-2)) 0] ))...
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
    q_plasma_cool       = P.gamma*e_charge* D2PM.T_L(Ntime)*mean(output.Gamma_n(Ntime,Nx-10:Nx)); D2PM.q_L(Ntime) = q_plasma_cool;
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
    xind = input.grid.i_xpoint(2);
    q_X = q_X_t(Ntime); 
    for Ntime = 1: length(output.time)
    D2PM.n_X(Ntime) = output.density(Ntime,xind);
    D2PM.n_L(Ntime) = mean(output.density(Ntime,Nx-30:Nx)); 
    D2PM.T_L(Ntime) = 1.5*output.temperature(Ntime,Nx)-0.5*output.temperature(Ntime,Nx-1);
    D2PM.T_X(Ntime) = output.temperature(Ntime,xind);
    D2PM.q_X(Ntime) = abs(q_X);
    D2PM.L   = L-L_core_SOL;
    % convective fraction
    f_cnv = (5*D2PM.n_X(Ntime)*e_charge*D2PM.T_X(Ntime)*output.velocity(Ntime,xind))/D2PM.q_X(Ntime);
    D2PM.f_cnv(Ntime) = max(min(0.9999,f_cnv),0);
    
    % power and momentum loss fractions with volume integral
    f_pwr = abs(-sum(delta_xcb(xind:end).*Source_Qtot(xind:end))/q_X);
    f_mom = abs(-sum(delta_xcb(xind:end).*(1./B_field(xind:end)*m.*density(Ntime,xind:end).*velocity(Ntime,xind:end).^2.*[diff(B_field_cb(xind:end))]./delta_xcb(xind:end) ...
                                    + Source_v(Ntime,xind:end)...
                                    + [input.numerics.viscosity*...
      (velocity(Ntime,xind+1:end)-2*velocity(Ntime,xind:end-1)+velocity(Ntime,xind-1:end-2)) 0] ))...
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
    q_plasma_cool       = P.gamma*e_charge* D2PM.T_L(Ntime)*mean(output.Gamma_n(Ntime,Nx-10:Nx)); D2PM.q_L(Ntime) = q_plasma_cool;
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
if P.plotoutput
 plot_div1d_output(output,input,proc)
end
if P.ploterror
 plot_div1d_error(output.X,error)
end
if P.plotprofiles
  plot_div1d_profiles(output, input,'avg',1,'plot',1,'figtight',1,'Ntime',P.Ntime,'print',1,...  );
     'save',0,'allpossiblefigures',0,'hold',1,'plot2PM',0,'xlimits',P.xlimits,'version',P.version,'FE',input.flux_expansion,'cratio',P.cratio)
end

%% output
out = struct;
out.proc = proc;
out.D2PM = D2PM;
out.D2PMwivi = D2PMwivi;
out.D2PMwovi = D2PMwovi;
out.Derror = error;

end

