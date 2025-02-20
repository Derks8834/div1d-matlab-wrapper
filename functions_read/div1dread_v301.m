function [output,input]= div1dread_v301(path)
% div1dread_v301(path) reads the output div1doutput.txt specified in path

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023


try % open the output file
    disp(['reading:  ',path])
    fid = fopen(path,'r');
catch
    input = struct; output = struct;
    disp('div1doutput.txt is not available');
    return
end
input = struct;
% read the numerics and physics parammeters
Text1 = textscan(fid, '%s', 1, 'delimiter','\n');
cell  = textscan(fid, '   Nx         = %u','delimiter','\n'); input.Nx    = cell2mat(cell(1));
cell  = textscan(fid, '   ntime      = %u','delimiter','\n'); input.ntime  = cell2mat(cell(1));
cell  = textscan(fid, '   nout       = %u','delimiter','\n'); input.nout  = cell2mat(cell(1));
cell  = textscan(fid, '   dxmin      = %f','delimiter','\n'); input.dxmin = cell2mat(cell(1));
cell  = textscan(fid, '   delta_t    = %f','delimiter','\n'); input.delta_t = cell2mat(cell(1));
cell  = textscan(fid, '   abstol     = %f','delimiter','\n'); input.abstol = cell2mat(cell(1));
cell  = textscan(fid, '   reltol     = %f','delimiter','\n'); input.reltol = cell2mat(cell(1));
cell  = textscan(fid, '   viscocity  = %f','delimiter','\n'); input.viscocity = cell2mat(cell(1));
cell  = textscan(fid, '   method     = %u','delimiter','\n'); input.method = cell2mat(cell(1));
cell  = textscan(fid, '   restart    = %s','delimiter','\n'); input.restart = cell{1}; %cell2mat(cell(1));
Text2 = textscan(fid, '%s', 1, 'delimiter','\n');
cell  = textscan(fid, '   gamma      = %f','delimiter','\n'); input.gamma  = cell2mat(cell(1));
cell  = textscan(fid, '   L          = %f','delimiter','\n'); input.L     = cell2mat(cell(1));
cell  = textscan(fid, '   sintheta   = %f','delimiter','\n'); input.sintheta= cell2mat(cell(1));
cell  = textscan(fid, '   mass       = %f','delimiter','\n'); input.mass     = cell2mat(cell(1));
cell  = textscan(fid, '   Gamma_X    = %f','delimiter','\n'); input.Gamma_X = cell2mat(cell(1));
cell  = textscan(fid, '   q_parX     = %f','delimiter','\n'); input.q_X   = cell2mat(cell(1));
cell  = textscan(fid, '   flux_exp   = %f','delimiter','\n'); input.flux_expansion = cell2mat(cell(1));
cell  = textscan(fid, '   initial_n  = %f','delimiter','\n'); input.initial_n = cell2mat(cell(1));
cell  = textscan(fid, '   initial_v  = %f','delimiter','\n'); input.initial_v = cell2mat(cell(1));
cell  = textscan(fid, '   initial_T  = %f','delimiter','\n'); input.initial_T = cell2mat(cell(1));
cell  = textscan(fid, '   initial_a  = %f','delimiter','\n'); input.initial_a = cell2mat(cell(1));
cell  = textscan(fid, '   density_ramp_rate       = %f','delimiter','\n'); input.density_ramp_rate = cell2mat(cell(1));
cell  = textscan(fid, '   energy_loss_ion         = %f','delimiter','\n'); input.energy_loss_ion = cell2mat(cell(1));
cell  = textscan(fid, '   neutral_residence_time  = %f','delimiter','\n'); input.neutral_residence_time = cell2mat(cell(1));
cell  = textscan(fid, '   redistributed_fraction  = %f','delimiter','\n'); input.redistributed_fraction = cell2mat(cell(1));
cell  = textscan(fid, '   recycling               = %f','delimiter','\n'); input.recycling = cell2mat(cell(1));
cell  = textscan(fid, '   num_impurities          = %f','delimiter','\n'); input.num_impurities = cell2mat(cell(1));
fmt = repmat(' %f ',1,input.num_impurities);
cell  = textscan(fid, strcat('   impurity_concentration  = ', fmt) ,input.num_impurities,'delimiter','\n'); input.impurity_concentration = cell2mat(cell(:));
cell  = textscan(fid, strcat('   impurity_Z              = ', fmt),input.num_impurities,'delimiter','\n'); input.impurity_Z = cell2mat(cell(:));
cell  = textscan(fid, '   gas_puff_source         = %f','delimiter','\n'); input.gas_puff_source = cell2mat(cell(1));
cell  = textscan(fid, '   gas_puff_location       = %f','delimiter','\n'); input.gas_puff_location = cell2mat(cell(1));
cell  = textscan(fid, '   gas_puff_width          = %f','delimiter','\n'); input.gas_puff_width = cell2mat(cell(1));
cell  = textscan(fid, '   elm_start_time          = %u','delimiter','\n'); input.elm_start_time = cell2mat(cell(1));
cell  = textscan(fid, '   elm_ramp_time           = %u','delimiter','\n'); input.elm_ramp_time = cell2mat(cell(1));
cell  = textscan(fid, '   elm_time_between        = %u','delimiter','\n'); input.elm_time_between = cell2mat(cell(1));
cell  = textscan(fid, '   elm_expelled_heat       = %f','delimiter','\n'); input.elm_expelled_heat = cell2mat(cell(1));
cell  = textscan(fid, '   elm_expelled_particles  = %f','delimiter','\n'); input.elm_expelled_particles = cell2mat(cell(1));
cell  = textscan(fid, '   switch_elm_series       = %u','delimiter','\n'); input.switch_elm_series = cell2mat(cell(1));
cell  = textscan(fid, '   gaussian_elm            = %f','delimiter','\n'); input.gaussian_elm = cell2mat(cell(1));
cell  = textscan(fid, '   radial_loss_factor      = %f','delimiter','\n'); input.radial_loss_factor = cell2mat(cell(1));
cell  = textscan(fid, '   radial_loss_gaussian    = %f','delimiter','\n'); input.radial_loss_gaussian = cell2mat(cell(1));
cell  = textscan(fid, '   radial_loss_width       = %f','delimiter','\n'); input.radial_loss_width = cell2mat(cell(1));
cell  = textscan(fid, '   radial_loss_location    = %f','delimiter','\n'); input.radial_loss_location = cell2mat(cell(1));
cell  = textscan(fid, '   switch_dyn_nu           = %u','delimiter','\n'); input.switch_dyn_nu = cell2mat(cell(1));
cell  = textscan(fid, '   switch_dyn_gas          = %u','delimiter','\n'); input.switch_dyn_gas = cell2mat(cell(1));
cell  = textscan(fid, '   switch_dyn_rec          = %u','delimiter','\n'); input.switch_dyn_rec = cell2mat(cell(1));
cell  = textscan(fid, '   switch_dyn_rad_los      = %u','delimiter','\n'); input.switch_dyn_rad_los = cell2mat(cell(1));
cell  = textscan(fid, '   switch_dyn_imp_con      = %u','delimiter','\n'); input.switch_car_con_prf = cell2mat(cell(1));
cell  = textscan(fid, '   switch_dyn_qpar         = %u','delimiter','\n'); input.switch_dyn_qpar = cell2mat(cell(1));
cell  = textscan(fid, '   switch_dyn_red_frc      = %u','delimiter','\n'); input.switch_dyn_red_frc = cell2mat(cell(1));

Nx = input.Nx;
% read the profile input data
% input.car_con_prf = zeros(Nx,1);
input.gas_puf_prf = zeros(Nx,1);
 % read the information line and skip
    Text = textscan(fid, '%s', 1, 'delimiter','\n');
    profile_inputs = textscan(fid,'%f %f %f ',Nx, 'delimiter','\n'); % 2 profiles now
    %Xdum = cell2mat(profile_inputs(:,1));
    input.x_grid = cell2mat(profile_inputs(:,1));
    input.gas_puf_prf = cell2mat(profile_inputs(:,2));
    input.flx_exp_prf = cell2mat(profile_inputs(:,3));
    input.B_field = input.flx_exp_prf;
    
Ntime = input.ntime/input.nout+1;
input.Ntime = Ntime;
% initialize data arrays
time = zeros(Ntime,1);
dyn_gas = zeros(Ntime,1);
dyn_nu = zeros(Ntime,1);
dyn_rec = zeros(Ntime,1);
dyn_rad_los = zeros(Ntime,1);
dyn_qpar = zeros(Ntime,1);
dyn_red_frc = zeros(Ntime,1);
dyn_imp_con = zeros(Ntime,input.num_impurities);

X = zeros(Nx,1);
density = zeros(Ntime,Nx);
velocity = zeros(Ntime,Nx);
temperature = zeros(Ntime,Nx);
neutral_density = zeros(Ntime,Nx);
Gamma_n = zeros(Ntime,Nx);
Gamma_mom = zeros(Ntime,Nx);
q_parallel = zeros(Ntime,Nx);
neutral_flux = zeros(Ntime,Nx);
Source_n = zeros(Ntime,Nx);
Source_v = zeros(Ntime,Nx);
Source_Q = zeros(Ntime,Nx);
Source_neutral = zeros(Ntime,Nx);

% get time data
for itime = 1: Ntime
    % read the first line containing the time
    cell = textscan(fid,'time        = %f','delimiter','\n'); time(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_gas     = %f','delimiter','\n'); dyn_gas(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_nu      = %f','delimiter','\n'); dyn_nu(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_rec     = %f','delimiter','\n'); dyn_rec(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_rad_los = %f','delimiter','\n'); dyn_rad_los(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_qparX   = %f','delimiter','\n'); dyn_qpar(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_red_frc = %f','delimiter','\n'); dyn_red_frc(itime) = cell2mat(cell(1));
    cell = textscan(fid,strcat('dyn_imp_con = ',fmt),'delimiter','\n'); dyn_imp_con(itime,1:input.num_impurities) = cell2mat(cell(:));
    
    % read the second line and skip
    Text = textscan(fid, '%s', 1, 'delimiter','\n');

    % now read the Nx data lines
    Output = textscan(fid,' %f %f %f %f %f %f %f %f %f %f %f %f %f ',Nx, 'delimiter','\n');
    X(:,1) = cell2mat(Output(:,1));
    density(itime,:) = cell2mat(Output(:,2));
    velocity(itime,:) = cell2mat(Output(:,3));
    temperature(itime,:) = cell2mat(Output(:,4));
    neutral_density(itime,:) = cell2mat(Output(:,5));
    Gamma_n(itime,:) = cell2mat(Output(:,6));
    Gamma_mom(itime,:) = cell2mat(Output(:,7));
    q_parallel(itime,:) = cell2mat(Output(:,8));
    neutral_flux(itime,:) = cell2mat(Output(:,9));
    Source_n(itime,:) = cell2mat(Output(:,10));
    Source_v(itime,:) = cell2mat(Output(:,11));
    Source_Q(itime,:) = cell2mat(Output(:,12));
    Source_neutral(itime,:) = cell2mat(Output(:,13));
end
cell = textscan(fid, '   cpu_time  = %f','delimiter','\n'); cpu_time = cell2mat(cell(1));
fclose(fid);

output = struct;
output.time         = time;
output.X            = X;
output.density      = density; 
output.velocity     = velocity;
output.temperature  = temperature;
output.neutral_density  = neutral_density;
output.Gamma_n      = Gamma_n;
output.Gamma_mom    = Gamma_mom;
output.q_parallel   = q_parallel;
output.neutral_flux = neutral_flux;
output.Source_n     = Source_n;
output.Source_v     = Source_v;
output.Source_Q     = Source_Q;
output.Source_neutral   = Source_neutral;
output.cpu_time = cpu_time;

% input = struct;
input.dyn_gas       = dyn_gas;
input.dyn_nu        = dyn_nu;
input.dyn_qpar      = dyn_qpar;
input.dyn_rad_los   = dyn_rad_los;
input.dyn_rec       = dyn_rec;
input.dyn_red_frc   = dyn_red_frc;
input.dyn_imp_con   = dyn_imp_con;
end