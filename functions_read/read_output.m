function [output,input] = read_output(path,varargin)
% [output,input]=read_output(dir,varargin)
% INPUT:        
%  dir      == string 'name of directory'    (e.g.== 'sID150683_v1.0.0_stat_maxqt_vis0.9_nu2_qu30_gv0_R74_fRd80/')
%  varargin == 'div1d_install_path','string' (def == '/home/emc/derks/Desktop')           
%           == 'dir','string'        (def == '/div1d/') (branch of div1dNOTE: should end with /)
% OUTPUT:
%   output = struct with entries:
%       [time, X, density, velocity, temperature, neutral_density, ...
%       Gamma_n, Gamma_mom, q_parallel, neutral_flux,...
%       Source_n, Source_v, Source_Q, Source_neutral, etc.]
%   input = struct with the inputs as USED by DIV1D. 

%% Handle parameters
if isunix
    D.div1d_install_path   = '/home/unix/derks/Desktop';
    D.dir                  = '';
else
    D.div1d_install_path   = '//Differ/Shares/Departments/Fusiefysica/IMM'; 
    D.dir                  = '';
end
D.runs_path            = '';
D.branch_path          ='/div1d/';
D.trymatfile           = true;  
P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
    if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end
% define path with output.txt
if ~isempty(P.dir)
path = [P.div1d_install_path,P.branch_path,'runs/',P.dir,'/div1d_output.txt'];
try pathc = strcat(path{:}); path = pathc; catch; end
end

if ~isempty(P.runs_path)
path = [ P.runs_path,'/div1d_output.txt'];
   disp(strcat('load: ', path))
try pathc = strcat(path{:}); path = pathc; catch; end
end

% try to open the equivalent .mat file
if P.trymatfile
    pathm = strrep(path,'.txt','.mat');
    try
        [output,input] =loadmatfile(pathm);
        disp('loaded .mat file');
        return
    catch
        disp('could not load .mat file');
    end
end
try
fid = fopen(path,'r');
try % read the code version that has written this output
    Text1 = textscan(fid,'git tag: %s','delimiter','\n');% git, tag: , v4.0.0 (version)
    version = Text1{1}; version = version{1};
    disp(['version: ',version]);
    P.version = version;
    versionfound = true;
    fclose(fid);
catch % otherwise try older/other versions
    %fclose(fid);
    disp('version not apparent, try all versions');
    allversions = {'v1.0.0','v2.0.0','v2.1.0','v3.0.0','v3.0.1','v3.0.2','v4.0.0','v4.0.1','v4.0.2','v6.0.0'};
    i = length(allversions);
    input = struct; output = struct;
    while isempty(fields(output))
        disp(['try:', allversions{i}])
        try
            [output, input] = read_version(path,allversions{i});
        catch
            disp('failed')
        end
        if i == 1; return; end
        i=i-1;
    end
    return
end
if versionfound
    [output, input] = read_version(path,P.version);
end

catch
disp('failed to open file, point at the .txt file');
[output, input] = read_version(path,'point');
end

return
end

function [output, input] = read_version(path,version)
% [output,input]=read_version(path,version) reads div1doutput
% e.g: [output,input]=read_version('/home/unix/derks/../div1d_output.txt','v4.0.0');
switch version
    case 'v1.0.0' % version of JP frankemoelle thesis
        [output,input]=div1dread_v100(path);
    case 'v2.0.0'  % version with flux expanion and dynamic inputs
        [output,input]=div1dread_v200(path);
    case 'v2.1.0'  % version with neutral background and correct recomb+friction
        [output,input]=div1dread_v210(path);
    case 'v3.0.0'  % version with multiple impurities (time-dependent)
        [output,input]=div1dread_v300(path);
    case 'v3.0.1'
        [output,input]=div1dread_v301(path);
    case 'v3.0.2'
        [output,input]=div1dread_v302(path);
    case 'v4.0.0'  % version now committed (full core-SOL + double X)
        [output,input]=div1dread_v400(path);
    case 'v4.0.1'  % version now committed (full core-SOL + double X)
        [output,input]=div1dread_v401(path);
    case 'v4.0.2'  % version now committed (full core-SOL + double X)
        [output,input]=div1dread_v402(path);
    case 'v6.0.0'  % version now committed (full core-SOL + double X)
        [output,input]=div1dread_v600(path);
    case 'point'

        [file,pathis] = uigetfile('*.txt');
        try
%         [output,input]=read_output(strcat(pathis,file));
        catch
        disp('file is not readable with known IO');
        end
    otherwise
        disp('this version is not available')
        input = struct;
        output = struct;
end
end

function [output,input]  = loadmatfile(path)
    load(path);
end
