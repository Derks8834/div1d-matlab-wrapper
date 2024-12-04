function [output,input] = div1dread_v401(path)
% div1dread_v401(path) reads the output div1doutput.txt specified in path

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
Text1 = textscan(fid,'git tag: %s','delimiter','\n');
version = Text1{1};
% we first read the namelists that are written
input = struct;
input.version = version;
input.numerics = read_namelist(fid, 'DIV1D_NUMERICS');
input.physics  = read_namelist(fid, 'DIV1D_PHYSICS' );
input.grid = read_namelist(fid, 'DIV1D_GRID'); 
input.physics.num_impurities = length(input.physics.impurity_concentration);

Nx = input.numerics.nx;   
Ntime = input.numerics.ntime/input.numerics.nout+1;
% add the legs for the core source profile
tmp = [zeros(1,input.grid.i_xpoint(1)-1) input.grid.core_source_profile zeros(1,Nx-input.grid.i_xpoint(2))] ; 
input.grid.core_source_profile = tmp;
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
dyn_q_core= zeros(Ntime,1);
dyn_gam_core = zeros(Ntime,1);
density = zeros(Ntime,Nx);
velocity = zeros(Ntime,Nx);
temperature = zeros(Ntime,Nx);
neutral_density = zeros(Ntime,Nx);
extern2sol_flux = zeros(Ntime,Nx);
Gamma_n = zeros(Ntime,Nx+1);
Gamma_mom = zeros(Ntime,Nx+1);
q_parallel = zeros(Ntime,Nx+1);
neutral_flux = zeros(Ntime,Nx+1);
Source_n = zeros(Ntime,Nx);
Source_v = zeros(Ntime,Nx);
Source_Q = zeros(Ntime,Nx);
Source_neutral = zeros(Ntime,Nx);

extern_neutral_flux = zeros(Ntime,3); 
sol2extern_flux= zeros(Ntime,5);
tar2extern_flux= zeros(Ntime,5);
extern_neutral_density= zeros(Ntime,5);
%end


 try
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
    cell = textscan(fid,'dyn_q_core  = %f','delimiter','\n'); dyn_q_core(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_gam_core= %f','delimiter','\n'); dyn_gam_core(itime) = cell2mat(cell(1));
    % read the second line with the headers and skip
    Text = textscan(fid, '%s', 1, 'delimiter','\n');

    % now read the Nx data lines for the variables and the sources at the cell centers
    Output = textscan(fid,'%f %f %f %f %f %f %f %f %f %f',Nx, 'delimiter','\n');
    X(:,1) = cell2mat(Output(:,1));
    density(itime,:) = cell2mat(Output(:,2));
    velocity(itime,:) = cell2mat(Output(:,3));
    temperature(itime,:) = cell2mat(Output(:,4));
    neutral_density(itime,:) = cell2mat(Output(:,5));
    extern2sol_flux(itime,:) = cell2mat(Output(:,6));
    Source_n(itime,:) = cell2mat(Output(:,7));
    Source_v(itime,:) = cell2mat(Output(:,8));
    Source_Q(itime,:) = cell2mat(Output(:,9));
    Source_neutral(itime,:) = cell2mat(Output(:,10));

    % read the next line with the headers and skip
    Text = textscan(fid, '%s', 1, 'delimiter','\n');
    % now read the Nx+1 data lines for the fluxes at the cell boundaries
    Output = textscan(fid,'%f %f %f %f %f ',Nx+1, 'delimiter','\n');
    Xcb(:,1) = cell2mat(Output(:,1));
    Gamma_n(itime,:) = cell2mat(Output(:,2));
    Gamma_mom(itime,:) = cell2mat(Output(:,3));
    q_parallel(itime,:) = cell2mat(Output(:,4));
    neutral_flux(itime,:) = cell2mat(Output(:,5));
    Text = textscan(fid,'%s',1,'delimiter','\n');
    Output = textscan(fid, '%f %f %f','delimiter','\n');
    extern_neutral_flux(itime,1:3) = cell2mat(Output(:)); 
    Text = textscan(fid,'%s',1,'delimiter','\n');
    Output = textscan(fid, '%f %f %f',5,'delimiter','\n');
    sol2extern_flux(itime,:) = cell2mat(Output(:,1));
    tar2extern_flux(itime,:) = cell2mat(Output(:,2));
    extern_neutral_density(itime,:) = cell2mat(Output(:,3));
end
cell = textscan(fid, '   cpu_time  = %f','delimiter','\n'); cpu_time = cell2mat(cell(1));
fclose(fid);

output = struct;
output.time         = time;
output.X            = input.grid.x;
output.Xcb          = input.grid.xcb; %input.grid.xcb;
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
output.extern2sol_flux = extern2sol_flux;

output.extern_neutral_flux  =  extern_neutral_flux;
output.sol2extern_flux      = sol2extern_flux;
output.tar2extern_flux      = tar2extern_flux;
output.extern_neutral_density= extern_neutral_density;

output.cpu_time = cpu_time;

catch
disp('could not read outputs')
output = struct;
end

input.dynamic.dyn_gas       = dyn_gas;
input.dynamic.dyn_nu        = dyn_nu;
input.dynamic.dyn_nb        = dyn_nb;
input.dynamic.dyn_qpar      = dyn_qpar;
input.dynamic.dyn_rad_los   = dyn_rad_los;
input.dynamic.dyn_rec       = dyn_rec;
input.dynamic.dyn_red_frc   = dyn_red_frc;
input.dynamic.dyn_imp_con   = dyn_imp_con;
input.dynamic.dyn_q_core    = dyn_q_core;
input.dynamic.dyn_gam_core  = dyn_gam_core;

try
pathm = strrep(path,'.txt','.mat');
save(pathm,"input","output",'-mat');
catch
disp('could not save .mat file');
end

end