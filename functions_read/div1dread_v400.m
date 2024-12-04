function [input,output] = div1dread_v400(path)
% div1dread_v400(path) reads the output div1doutput.txt specified in path

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

% read the code version that has written this output
Text1 = textscan(fid,'git tag: %s','delimiter','\n');% git, tag: , v4.0.0 (version)
version = Text1{1};
% we first read the namelists that are written
input = struct;
input.numerics = read_namelist(fid, 'DIV1D_NUMERICS');
input.physics  = read_namelist(fid, 'DIV1D_PHYSICS' );
%input.internal = read_namelist(fid, 'DIV1D_INTERNAL'); 

Nx = input.numerics.nx;
% read the profile input data
input.gas_puf_prf = zeros(Nx,1);
% read the information line and skip
Text = textscan(fid, '%s', 1, 'delimiter','\n'); % X [m] gas_puff_prf [] B_field [fraction]
profile_inputs = textscan(fid,'%f %f %f ',Nx, 'delimiter','\n'); % 2 profiles now
input.x_grid = cell2mat(profile_inputs(:,1));
input.gas_puf_prf = cell2mat(profile_inputs(:,2));
input.B_field = cell2mat(profile_inputs(:,3));
    
Ntime = input.numerics.ntime/input.numerics.nout+1;

% initialize data arrays
time = zeros(Ntime,1);
dyn_gas = zeros(Ntime,1);
dyn_nu = zeros(Ntime,1);
dyn_nb = zeros(Ntime,1);
dyn_rec = zeros(Ntime,1);
dyn_rad_los = zeros(Ntime,1);
dyn_qpar = zeros(Ntime,1);
dyn_red_frc = zeros(Ntime,1);
dyn_imp_con = zeros(Ntime,max(input.physics.num_impurities,1));

X = zeros(Nx,1);
Xcb = zeros(Nx+1,1); % note that the cell boundary count is 1: Nx+1 instead of 0:NX as in the fortran code
density = zeros(Ntime,Nx);
velocity = zeros(Ntime,Nx);
temperature = zeros(Ntime,Nx);
neutral_density = zeros(Ntime,Nx);
Gamma_n = zeros(Ntime,Nx+1);
Gamma_mom = zeros(Ntime,Nx+1);
q_parallel = zeros(Ntime,Nx+1);
neutral_flux = zeros(Ntime,Nx+1);
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
    cell = textscan(fid,'dyn_nb      = %f','delimiter','\n'); dyn_nb(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_rec     = %f','delimiter','\n'); dyn_rec(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_rad_los = %f','delimiter','\n'); dyn_rad_los(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_qparX   = %f','delimiter','\n'); dyn_qpar(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_red_frc = %f','delimiter','\n'); dyn_red_frc(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_imp_con = %f %f %f %f %f','delimiter','\n'); dyn_imp_con(itime,1:5) = cell2mat(cell(:));
    
    % read the second line with the headers and skip
    Text = textscan(fid, '%s', 1, 'delimiter','\n');

    % now read the Nx data lines for the variables and the sources at the cell centers
    Output = textscan(fid,' %f %f %f %f %f %f %f %f %f ',Nx, 'delimiter','\n');
    X(:,1) = cell2mat(Output(:,1));
    density(itime,:) = cell2mat(Output(:,2));
    velocity(itime,:) = cell2mat(Output(:,3));
    temperature(itime,:) = cell2mat(Output(:,4));
    neutral_density(itime,:) = cell2mat(Output(:,5));
    Source_n(itime,:) = cell2mat(Output(:,6));
    Source_v(itime,:) = cell2mat(Output(:,7));
    Source_Q(itime,:) = cell2mat(Output(:,8));
    Source_neutral(itime,:) = cell2mat(Output(:,9));

    % read the next line with the headers and skip
    Text = textscan(fid, '%s', 1, 'delimiter','\n');
    % now read the Nx+1 data lines for the fluxes at the cell boundaries
    Output = textscan(fid,' %f %f %f %f %f ',Nx+1, 'delimiter','\n');
    Xcb(:,1) = cell2mat(Output(:,1));
    Gamma_n(itime,:) = cell2mat(Output(:,2));
    Gamma_mom(itime,:) = cell2mat(Output(:,3));
    q_parallel(itime,:) = cell2mat(Output(:,4));
    neutral_flux(itime,:) = cell2mat(Output(:,5));

end
cell = textscan(fid, '   cpu_time  = %f','delimiter','\n'); cpu_time = cell2mat(cell(1));
fclose(fid);

output = struct;
output.time         = time;
output.X            = X;
output.Xcb          = Xcb;
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

input.dyn_gas       = dyn_gas;
input.dyn_nu        = dyn_nu;
input.dyn_nb        = dyn_nb;
input.dyn_qpar      = dyn_qpar;
input.dyn_rad_los   = dyn_rad_los;
input.dyn_rec       = dyn_rec;
input.dyn_red_frc   = dyn_red_frc;
input.dyn_imp_con   = dyn_imp_con;

end