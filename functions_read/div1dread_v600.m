function [output,input] = div1dread_v600(path)
%[output,input]= div1dread_v600(path) reads the output div1doutput.txt specified in path

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% Aug 2024

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
%input.grid.vesrz = reshape(input.grid.vesrz,3,length(input.grid.vesrz)/3);
%         indiv{i_set}.geom.vslr = geometry.vesrz(:,1)';
%         indiv{i_set}.geom.vslz = geometry.vesrz(:,2)';
%         indiv{i_set}.geom.nr = geometry.normal_vector(1,:);
%         indiv{i_set}.geom.nz = geometry.normal_vector(2,:);
% below should be generated from read files
%         indiv{i_set}.geom.ir = geometry.ir(1,:);
%         indiv{i_set}.geom.iz = geometry.iz(1,:);
%         indiv{i_set}.geom.or = geometry.or(1,:);
%         indiv{i_set}.geom.oz = geometry.oz(1,:);
%         indiv{i_set}.geom.romrz = romrz;



Nx = input.numerics.nx;   
Ntime = ceil((input.numerics.ntime-1)/input.numerics.nout)+1;
input.grid.prf_imp_con = reshape(input.grid.prf_imp_con,Nx,input.physics.num_impurities);

% add the legs for the core source profile
% if input.physics.l_core_sol > 0
% tmp = [zeros(1,input.grid.i_xpoint(1)-1) input.grid.core_source_profile zeros(1,Nx-input.grid.i_xpoint(2))] ; 
% input.grid.core_source_profile = tmp;
% else 
%     input.grid.core_source_profile = nan;
% end
% initialize data arrays
time = zeros(Ntime,1);
%dyn = struct;
dyn_gas = zeros(Ntime,1);
dyn_nu = zeros(Ntime,1);
dyn_nb = zeros(Ntime,5);
dyn_mb = zeros(Ntime,5);
dyn_rec = zeros(Ntime,1);
dyn_qpar = zeros(Ntime,1);
dyn_imp_con = zeros(Ntime,max(input.physics.num_impurities,1));
dyn_q_core= zeros(Ntime,1);
dyn_gam_core = zeros(Ntime,1);
dyn_neutral_puff = zeros(Ntime,5); 
dyn_molecule_puff = zeros(Ntime,5);
dyn_core_fuelling = zeros(Ntime,1);
dyn_core_neutral_density = zeros(Ntime,1);

density = zeros(Ntime,Nx);
velocity = zeros(Ntime,Nx);
temperature = zeros(Ntime,Nx);
neutral_density = zeros(Ntime,Nx);
neutral_velocity = zeros(Ntime,Nx);
molecule  = zeros(Ntime,Nx);
extern2sol_flux = zeros(Ntime,Nx);
extern2sol_mol = zeros(Ntime,Nx);
core2sol_flux = zeros(Ntime,Nx);
core2sol_mol = zeros(Ntime,Nx);
sol2extern_ion_flux = zeros(Ntime,Nx);
Gamma_n = zeros(Ntime,Nx+1);
Gamma_mom = zeros(Ntime,Nx+1);
q_parallel = zeros(Ntime,Nx+1);
Gamma_neutral = zeros(Ntime,Nx+1);
Gamma_mom_neutral = zeros(Ntime,Nx+1);
Gamma_molecule =  zeros(Ntime,Nx+1);

Source_n = zeros(Ntime,Nx);
Source_v = zeros(Ntime,Nx);
Source_Q = zeros(Ntime,Nx);
Source_neutral = zeros(Ntime,Nx);
Source_vn = zeros(Ntime,Nx);
Source_molecule = zeros(Ntime,Nx);

% external chambers
extern_neutral_flux = zeros(Ntime,3); 
sol2extern_flux= zeros(Ntime,5);
tar2extern_flux= zeros(Ntime,5);
extern2core_flux = zeros(Ntime,5);
extern_neutral_density= zeros(Ntime,5);
extern_molecule_flux = zeros(Ntime,3); 
sol2extern_mol= zeros(Ntime,5);
tar2extern_mol = zeros(Ntime,5);
extern2core_mol = zeros(Ntime,5);
extern_molecule_density= zeros(Ntime,5);

% core
core_density = zeros(Ntime,1); 
Source_core = zeros(Ntime,1); 
core_neutral_density = zeros(Ntime,1); 
Gamma_core2sol = zeros(Ntime,1); 
sol2core_flux = zeros(Ntime,1); 
sol2core_mol = zeros(Ntime,1); 
sum_extern2core_flux = zeros(Ntime,1);
sum_extern2core_mol = zeros(Ntime,1);

% strange ducks
sum_sol2extern_ion_flux = zeros(Ntime,1); 
sum_sol2extern_ion_mol = zeros(Ntime,1); 
  fs18 = fstr(18);
%  try
% get time data
try 
for itime = 1: Ntime
    % read the first line containing the time dependent inputs
    cell = textscan(fid,'time        = %f','delimiter','\n'); time(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_gas     = %f','delimiter','\n'); dyn_gas(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_nu      = %f','delimiter','\n'); dyn_nu(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_nb      = %f %f %f %f %f','delimiter','\n'); dyn_nb(itime,1:5) = cell2mat(cell(:));
    cell = textscan(fid,'dyn_mb      = %f %f %f %f %f','delimiter','\n'); dyn_mb(itime,1:5) = cell2mat(cell(:));
    cell = textscan(fid,'dyn_rec     = %f','delimiter','\n'); dyn_rec(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_qparX   = %f','delimiter','\n'); dyn_qpar(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_imp_con = %f %f %f %f %f','delimiter','\n'); dyn_imp_con(itime,1:5) = cell2mat(cell(:));
    cell = textscan(fid,'dyn_q_core  = %f','delimiter','\n'); dyn_q_core(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_gam_core= %f','delimiter','\n'); dyn_gam_core(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_neutral_puff = %f %f %f %f %f','delimiter','\n'); dyn_neutral_puff(itime,1:5) = cell2mat(cell(:));
    cell = textscan(fid,'dyn_molecule_puff = %f %f %f %f %f','delimiter','\n'); dyn_molecule_puff(itime,1:5) = cell2mat(cell(:));
    cell = textscan(fid,'dyn_core_fuelling = %f','delimiter','\n'); dyn_core_fuelling(itime) = cell2mat(cell(1));
    cell = textscan(fid,'dyn_core_neutral_density = %f','delimiter','\n'); dyn_core_neutral_density(itime) = cell2mat(cell(1));


    % read the header for SOL quantities sources and fluxes and skip 
    Text = textscan(fid, '%s', 1, 'delimiter','\n');

    % now read the Nx data lines for the variables and the sources at the cell center
    Output = textscan(fid,fs18,Nx, 'delimiter','\n');
    it = 1;
    X(:,1) = cell2mat(Output(:,it));         it = it +1;
    density(itime,:) = cell2mat(Output(:,it));it = it +1;
    velocity(itime,:) = cell2mat(Output(:,it));it = it +1;
    temperature(itime,:) = cell2mat(Output(:,it));it = it +1;
    neutral_density(itime,:) = cell2mat(Output(:,it));it = it +1;
    neutral_velocity(itime,:) = cell2mat(Output(:,it));it = it +1;
    molecule(itime,:) = cell2mat(Output(:,it));it = it +1;

    extern2sol_flux(itime,:) = cell2mat(Output(:,it));it = it +1;
    extern2sol_mol(itime,:) = cell2mat(Output(:,it));it = it +1;
    core2sol_flux(itime,:) = cell2mat(Output(:,it));it = it +1;
    core2sol_mol(itime,:) = cell2mat(Output(:,it));it = it +1;
    sol2extern_ion_flux(itime,:) = cell2mat(Output(:,it));it = it +1;
  
    Source_n(itime,:) = cell2mat(Output(:,it));it = it +1;
    Source_v(itime,:) = cell2mat(Output(:,it));it = it +1;
    Source_Q(itime,:) = cell2mat(Output(:,it));it = it +1;
    Source_neutral(itime,:) = cell2mat(Output(:,it));it = it +1;
    Source_vn(itime,:) = cell2mat(Output(:,it));it = it +1;
    Source_molecule(itime,:) = cell2mat(Output(:,it));%it = it +1;
%(x(i), density(i), velocity(i), temperature(i), neutral(i), neutral_velocity(i), molecule(i), & ! 7!
%				extern2sol_flux(i), extern2sol_mol(i), core2sol_flux(i), core2sol_mol(i), sol2extern_ion_flux(i),  & ! 5
%				 Source_n(i), Source_v(i), Source_Q(i), Source_neutral(i), Source_vn(i), Source_molecule(i), i=1,Nx ) ! 6 -> 18
    % read the next line with the headers and skip fluxes
    Text = textscan(fid, '%s', 1, 'delimiter','\n');
    % now read the Nx+1 data lines for the fluxes at the cell boundaries
    Output = textscan(fid,'%f %f %f %f %f %f %f ',Nx+1, 'delimiter','\n');
    Xcb(:,1) = cell2mat(Output(:,1));
    Gamma_n(itime,:) = cell2mat(Output(:,2));
    Gamma_mom(itime,:) = cell2mat(Output(:,3));
    q_parallel(itime,:) = cell2mat(Output(:,4));
    Gamma_neutral(itime,:) = cell2mat(Output(:,5));
    Gamma_mom_neutral(itime,:) = cell2mat(Output(:,6));
    Gamma_molecule(itime,:) = cell2mat(Output(:,7));
    
    % extern neutrals all in units [D]
    Text = textscan(fid,'%s',1,'delimiter','\n'); % extern neutral flux
    Output = textscan(fid, '%f %f %f','delimiter','\n');
    extern_neutral_flux(itime,1:3) = cell2mat(Output(:)); 
    Text = textscan(fid,'%s',1,'delimiter','\n');
    Output = textscan(fid, '%f %f %f %f',5,'delimiter','\n');
    sol2extern_flux(itime,:) = cell2mat(Output(:,1));
    tar2extern_flux(itime,:) = cell2mat(Output(:,2));
    extern2core_flux(itime,:) = cell2mat(Output(:,3));
    extern_neutral_density(itime,:) = cell2mat(Output(:,4));

    % extern molecules all in units [D2]
    Text = textscan(fid,'%s',1,'delimiter','\n'); % extern molecule flux
    Output = textscan(fid, '%f %f %f','delimiter','\n');
    extern_molecule_flux(itime,1:3) = cell2mat(Output(:)); 
    Text = textscan(fid,'%s',1,'delimiter','\n'); % sol2ex tar2ex ext dens
    Output = textscan(fid, '%f %f %f %f',5,'delimiter','\n');
    sol2extern_mol(itime,:) = cell2mat(Output(:,1))';
    tar2extern_mol(itime,:) = cell2mat(Output(:,2))';
    extern2core_mol(itime,:) = cell2mat(Output(:,3))';
    extern_molecule_density(itime,:) = cell2mat(Output(:,4));

    % the ducks
    Text = textscan(fid,'%s',1,'delimiter','\n'); % sum sol2extern ion flux and mol
    Output = textscan(fid, '%f %f',1,'delimiter','\n');
    sum_sol2extern_ion_flux(itime,:) = cell2mat(Output(1,1));
    sum_sol2extern_ion_mol(itime,:) = cell2mat(Output(1,2));

    % core  all in units [D]
    Text = textscan(fid,'%s',1,'delimiter','\n'); % density, source, neutral, gamma, sol2core
    Output = textscan(fid, '%f %f %f %f %f %f %f %f',1,'delimiter','\n');
    core_density(itime,:) = cell2mat(Output(1,1));
    Source_core(itime,:) = cell2mat(Output(1,2));
    core_neutral_density(itime,:) = cell2mat(Output(1,3));
    Gamma_core2sol(itime,:) = cell2mat(Output(1,4));
    sol2core_flux(itime,:) =  cell2mat(Output(1,5));
    sol2core_mol(itime,:) = cell2mat(Output(1,6));
    sum_extern2core_flux(itime,:) = cell2mat(Output(1,7));
    sum_extern2core_mol(itime,:) = cell2mat(Output(1,8));


    % start over

end
cell = textscan(fid, '   cpu_time  = %f','delimiter','\n'); cpu_time = cell2mat(cell(1));
catch
disp('did not read all timesteps')
end
fclose(fid);

output = struct;
output.time         = time;
% sol
output.X            = input.grid.x;
output.Xcb          = input.grid.xcb; %input.grid.xcb;
output.density      = density; 
output.velocity     = velocity;
output.neutral_velocity = neutral_velocity;
output.temperature  = temperature;
output.neutral_density  = neutral_density;
output.molecule = molecule;

output.Gamma_n      = Gamma_n;
output.Gamma_mom    = Gamma_mom;
output.q_parallel   = q_parallel;
output.Gamma_neutral = Gamma_neutral;
output.Gamma_mom_neutral = Gamma_mom_neutral;
output.Gamma_molecule = Gamma_molecule;

output.Source_n     = Source_n;
output.Source_v     = Source_v;
output.Source_Q     = Source_Q;
output.Source_neutral  = Source_neutral;
output.Source_vn       = Source_vn;
output.Source_molecule = Source_molecule;
% neutral fluxes from cross-channel transport
output.extern2sol_flux      = extern2sol_flux; 
output.extern2sol_mol = extern2sol_mol;
output.core2sol_flux = core2sol_flux;
output.core2sol_mol = core2sol_mol;

% atoms [D]
output.extern_neutral_flux  = extern_neutral_flux;
output.sol2extern_flux      = sol2extern_flux;
output.tar2extern_flux      = tar2extern_flux;
output.extern2core_flux     = extern2core_flux;
output.extern_neutral_density  = extern_neutral_density;
output.sum_sol2extern_ion_flux = sum_sol2extern_ion_flux;
% molecules [D2]
output.extern_molecule_flux  = extern_molecule_flux;
output.sol2extern_mol        = sol2extern_mol;
output.tar2extern_mol        = tar2extern_mol;
output.extern2core_mol       = extern2core_mol;
output.extern_molecule_density = extern_molecule_density;
output.sum_sol2extern_ion_mol  = sum_sol2extern_ion_mol;

% core [D]
output.core_density = core_density;
output.Source_core = Source_core;
output.core_neutral_density = core_neutral_density;
output.Gamma_core2sol = Gamma_core2sol;
output.sol2core_flux = sol2core_flux;
output.sol2core_mol  = sol2core_mol;
output.sum_extern2core_flux = sum_extern2core_flux;
output.sum_extern2core_mol  = sum_extern2core_mol;

output.cpu_time = cpu_time;

% catch
% disp('could not read outputs')
% output = struct;
% end

input.dynamic.dyn_gas       = dyn_gas;
input.dynamic.dyn_nu        = dyn_nu;
input.dynamic.dyn_nb        = dyn_nb;
input.dynamic.dyn_mb        = dyn_mb;
input.dynamic.dyn_qpar      = dyn_qpar;
input.dynamic.dyn_rec       = dyn_rec;
input.dynamic.dyn_imp_con   = dyn_imp_con;
input.dynamic.dyn_q_core    = dyn_q_core;
input.dynamic.dyn_gam_core  = dyn_gam_core;
input.dynamic.dyn_neutral_puff      = dyn_neutral_puff;
input.dynamic.dyn_molecule_puff     = dyn_molecule_puff;
input.dynamic.dyn_core_fuelling     = dyn_core_fuelling;
input.dynamic.dyn_core_neutral_density = dyn_core_neutral_density;

try
pathm = strrep(path,'.txt','.mat');
save(pathm,"input","output",'-mat');
catch
disp('could not save .mat file');
end

    function [out] = fstr(N)
        % input N number of %f that have to be appended
        % out = '%f %f %f ' N times

        out = '%f';
        if N > 1
        for i = 1:N-1
        out = strcat(out,' %f') ;
        end
        end
    end

end