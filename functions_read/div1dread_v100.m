function  [output,input] = div1dread_v100(path)
%div1dread_v100(path) reads the output div1doutput.txt specified in path

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
% read the numerics and physics parammeters
Text1 = textscan(fid, '%s', 1, 'delimiter','\n');
cell  = textscan(fid, '   Nx         = %u','delimiter','\n');
Nx    = cell2mat(cell(1));
cell  = textscan(fid, '   ntime      = %u','delimiter','\n');
ntime = cell2mat(cell(1));
cell  = textscan(fid, '   nout       = %u','delimiter','\n');
nout  = cell2mat(cell(1));
cell  = textscan(fid, '   dxmin      = %f','delimiter','\n');
dxmin = cell2mat(cell(1));
Text2 = textscan(fid, '%s', 1, 'delimiter','\n');
cell  = textscan(fid, '   L          = %f','delimiter','\n');
L     = cell2mat(cell(1));
cell   = textscan(fid, '   q_parX     = %f','delimiter','\n');
q_X   = cell2mat(cell(1));
cell  = textscan(fid, '   initial_n  = %f','delimiter','\n');
initial_n = cell2mat(cell(1));
Ntime = ntime/nout+1;

input = struct;
input.Nx = Nx;
input.ntime = ntime;
input.nout = nout;
input.dxmin = dxmin;
input.L = L;
input.q_X =q_X;
input.initial_n = initial_n;
input.Ntime = Ntime;

% initialize data arrays
time = zeros(Ntime,1);
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
    TimeOut = textscan(fid, ' time =  %f','delimiter','\n');
    time(itime) = cell2mat(TimeOut(1));
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
end