function [flag] =write_input(numpar,phypar,varargin)
% write_input(numpar,phypar,varargin)
% - phypar and numpar are structs and turned into namelists for DIV1D
% - varargin: time- and space-dependent inputs are parsed via ".dat" files.
% varargin arguments in dat struct
%   'nu'        ,vec (with vec[1:ntime])
%   'R'         ,vec "
%   'gas'       ,vec "
%   'RL'        ,vec "
%   'fr'        ,vec "
%   'qpar'      ,vec "
%   'imp'       ,vec (with [5,ntime]) 
%   'car_frc'   ,vec (with [5,1:Nx])
%  
% NOTE: run_program.m must be instructed to use/copy ".dat' files with <varargin 'datinput',1> 
%
% Output: flag ( = true if everything worked fine)
%
% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

D = struct;
D.dat = struct;
%% Handle Inputs
P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end

try
system('bash -c "rm ./input.txt"');
catch
end
try
system('bash -c "rm ./*.dat"');
catch
end

%% output flag
flag = true;

%% write namelists
fid=fopen('input.txt','w');

write_fortran_namelist(fid,numpar,'div1d_numerics');
write_fortran_namelist(fid,phypar,'div1d_physics');
fclose(fid);


%%  Write time and space dependent inputs.
 ntime   = numpar.ntime; % number of time steps in input
 Nx      = numpar.Nx; % number of cells along leg
 dat=P.dat;

%Write .dat files
     try 
    if min(phypar.impurity_concentration) < -0.5 %'impurities';
           dyn_imp_con = dat.impurity_concentration; 
           
        if length(dyn_imp_con) >1
            if length(dyn_imp_con) == ntime
            filename = 'dyn_imp_con.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            strtmp = repmat('%5.10f ', 1,5); % ntime
            fprintf(fid,strcat(strtmp,' \n'),dyn_imp_con);
            fclose(fid);
            else
                disp('warning: dyn_imp_con.dat and ntime have a discrepancy')
                flag = false;
            end
        end
    end
     catch end
     try
    if phypar.initial_n < 0 %'upstream density';
           dyn_nu = dat.initial_n; 
        if length(dyn_nu) >1
            if length(dyn_nu) == ntime
            filename = 'dyn_nu.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',dyn_nu);
            fclose(fid);
            else
                disp('warning: dyn_nu.dat and ntime have a discrepancy')
                flag = false;
            end
        end
    end
     catch end
     try
    if phypar.initial_a < 0 %'neutral background';
                dyn_nb = dat.initial_a; 
        if length(dyn_nb) >1
            if length(dyn_nb) == ntime
            filename = 'dyn_nb.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',dyn_nb);
            fclose(fid);
            else
                disp('warning: dyn_nb.dat and ntime have a discrepancy')
                flag = false;
            end
        end
    end
     catch end
     try
    if phypar.q_par_X < 0 %'qpar'; 
                dyn_qpar = dat.q_parX_; 
        if length(dyn_qpar) >1
            if length(dyn_qpar) == ntime
            filename = 'dyn_qpar.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',dyn_qpar);
            fclose(fid);
            else
                disp('warning: dyn_qpar.dat and ntime have a discrepancy')
                flag = false;
            end
        end
    end
     catch
     end
    try
    if phypar.redistributed_fraction < 0 %'redistribution fraction';
                dyn_red_frc = dat.redistributed_fraction; 
        if length(dyn_red_frc) >1
            if length(dyn_red_frc) == ntime
            filename = 'dyn_red_frc.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',dyn_red_frc);
            fclose(fid);
            else
                disp('warning: dyn_red_frc.dat and ntime have a discrepancy')
                flag = false;
            end
        end
    end
    catch end
    try
    if phypar.recycling < 0 %'recycling';
             dyn_rec = dat.recycling;
        if length(dyn_rec) >1
             if length(dyn_rec) == ntime
            filename = 'dyn_rec.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',dyn_rec);
            fclose(fid);
           else
                disp('warning: dyn_rec.dat and ntime have a discrepancy')
                flag = false;
            end
        end
    end
    catch end
    try
    if phypar.gas_puff_source < 0 %'gas puff source';
                dyn_gas = dat.gas_puff_source;
        if length(dyn_gas) >1
             if length(dyn_gas) == ntime+2
            filename = 'dyn_gas.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',dyn_gas);
            fclose(fid);
           else
                disp('warning: dyn_gas.dat and ntime have a discrepancy')
                flag = false;
            end
        end
    end
    catch end
  

  if phypar.sintheta < 0 % sintheta ;
                dyn_sintheta = dat.sintheta;
        if length(dyn_sintheta) >1
             if length(dyn_sintheta) == Nx
            filename = 'sintheta.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',dyn_sintheta);
            fclose(fid);
           else
                disp('warning: sintheta.dat and Nx have a discrepancy')
                flag = false;
            end
        end
  end
 if phypar.major_radius < 0 % width outer midplane;
                prf_major_radius = dat.major_radius;
        if length(prf_major_radius) >1
            if length(prf_major_radius) == Nx
            filename = 'major_radius.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',prf_major_radius);
            fclose(fid);
           else
                disp('warning: major_radius.dat and Nx have a discrepancy')
                flag = false;
            end
        end
        try
        prf_major_height = dat.major_height;
        if length(prf_major_height) >1
            if length(prf_major_height) == Nx
            filename = 'major_height.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',prf_major_height);
            fclose(fid);
           else
                disp('warning: major_height.dat and Nx have a discrepancy')
                flag = false;
            end
        end
        catch
        disp( 'failed to write Z to DIV1D for plotting purposes')
        end
        try
        prf_sol_normal = dat.sol_normal;
         if length(prf_sol_normal) >1
            if length(prf_sol_normal) == Nx
            filename = 'sol_normal.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f  %5.10f \n',prf_sol_normal);
            fclose(fid);
           else
                disp('warning: sol_normal.dat and Nx have a discrepancy')
            end
        end
        catch
        disp( 'failed to write sol_normal.dat to DIV1D for plotting purposes')
        end
        try
         [n,m]  = size(dat.vesrz);
         msum = sum(dat.vesrz(:,m));
         if m ==3 && msum ==9 % 9 is number of markers on wall
        vessel = dat.vesrz;
         if length(vessel) >1
            if length(vessel) < Nx
            filename = 'vessel_r.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%d \n',n);
            for i_n = 1:n
                fprintf(fid,'%5.10f \n',vessel(i_n,1));                
            end
            fclose(fid);

            filename = 'vessel_z.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%d \n',n);
            for i_n = 1:n
                fprintf(fid,'%5.10f \n',vessel(i_n,2));                
            end
            fclose(fid);
            
            filename = 'vessel_c.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%d \n',n);
            for i_n = 1:n
                fprintf(fid,'%5.10f \n',vessel(i_n,3));                
            end
            fclose(fid);

            else
                disp('warning: vessel_rzc.dat is larger than Nx, the max..')
            end
         end
         else
             disp('vesrz should contain 3 columns with r value, z value and a 0/1 indicating a cut iff 1, for a total of 9 cutts')
         end
        catch
        disp( 'failed to write sol_normal.dat to DIV1D for plotting purposes')
        end
      
 end
 if phypar.flux_expansion < 0 % get array
        prf_B_field = dat.B_field;
        if length(prf_B_field) >1
             if length(prf_B_field) == Nx
            filename = 'B_field.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',prf_B_field);
            fclose(fid);
           else
                disp('warning: B_field.dat and Nx have a discrepancy')
                flag = false;
            end
        end
 end
  if phypar.trans_expansion < 0 % get array
        prf_B_trans = dat.B_trans;
        if length(prf_B_trans) >1
             if length(prf_B_trans) == Nx
            filename = 'B_trans.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',prf_B_trans);
            fclose(fid);
           else
                disp('warning: B_trans.dat and Nx have a discrepancy')
                flag = false;
            end
        end
 end
  if phypar.Gamma_core< 0 % width outer midplane;
                dyn_gamma_core = dat.gamma_core;
        if length(dyn_gamma_core) >1
            if length(dyn_gamma_core) == ntime+1    
            filename = 'dyn_gamma_core.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            fprintf(fid,'%5.10f \n',dyn_gamma_core);
            fclose(fid);
           else
                disp('warning: dyn_gamma_core.dat and ntime have a discrepancy')
                flag = false;
            end
        end
  end

   if min(phypar.initial_nb)< -0.5  % neutral background densities;
                dyn_nb = dat.initial_nb;
        if length(dyn_nb) >1
            logtmp = (size(dyn_nb) == [5,ntime+1]);
             if logtmp(1) &&  logtmp(2)  % length(dyn_nb) == ntime+1    
            filename = 'dyn_nb.dat'; permissions = 'w';
            fid = fopen(filename,permissions);
            strtmp = repmat('%5.10f ', 1,5); % ntime
            fprintf(fid,strcat(strtmp,' \n'),dyn_nb);
            fclose(fid);
           else
                disp('warning: dyn_nb.dat and ntime have a discrepancy')
                flag = false;
            end
        end
   end

   if min(phypar.initial_mb)< -0.5  % neutral background densities;
            dyn_mb = dat.initial_mb;
    if length(dyn_mb) >1
        logtmp = (size(dyn_mb) == [5,ntime+1]);
        if logtmp(1) &&  logtmp(2)  %length(dyn_mb) == ntime+1    
        filename = 'dyn_mb.dat'; permissions = 'w';
        fid = fopen(filename,permissions);
        strtmp = repmat('%5.10f ', 1,5); % ntime
        fprintf(fid,strcat(strtmp,' \n'),dyn_mb);
        fclose(fid);
       else
            disp('warning: dyn_mb.dat and ntime have a discrepancy')
            flag = false;
        end
    end
   end
  if isfield(phypar,'puff_rate_molecule');
   if phypar.puff_rate_molecule(1) == -1 %-0.5  % molecule puffing;
            dyn_molecule_puff = dat.molecule_puff;
    if length(dyn_molecule_puff) >1
        logtmp = (size(dyn_molecule_puff) == [5,ntime+1]);
        if logtmp(1) && logtmp(2) %length(dyn_molecule_puff) == ntime+1    
        filename = 'molecule_puff.dat'; permissions = 'w';
        fid = fopen(filename,permissions);
        strtmp = repmat('%5.10f ', 1,5); % ntime
        fprintf(fid,strcat(strtmp,' \n'),dyn_molecule_puff);
        fclose(fid);
       else
            disp('warning: molecule_puff.dat and ntime have a discrepancy')
            flag = false;
        end
    end

   end
  end
    if isfield(phypar,'neutral_puff');
    if phypar.puff_rate_neutral(1) == -1 % -0.5  % neutral atom puffing;
            dyn_neutral_puff = dat.neutral_puff;
    if length(dyn_neutral_puff) >1
        logtmp = (size(dyn_neutral_puff) == [5,ntime+1]);
        if logtmp(1) &&  logtmp(2)  %length(dyn_molecule_puff) == ntime+1    
        filename = 'neutral_puff.dat'; permissions = 'w';
        fid = fopen(filename,permissions);
        strtmp = repmat('%5.10f ', 1,5); % ntime
        fprintf(fid,strcat(strtmp,' \n'),dyn_neutral_puff);
        fclose(fid);
       else
            disp('warning: neutral_puff.dat and ntime have a discrepancy')
            flag = false;
        end
    end
   end
    end
  
end
