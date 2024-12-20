function [geometry] = plot_div1d_polgeom(indiv, varargin)
%PLOT_DIV1D_POLGEOM Summary of this function goes here
% function [h] = plot_div1d_polgeom(indiv, varargin)
% plots the poloidal geometry of div1d
% input: indiv (input struct div1d)
%        varargin ('name',value) pairs with default values:
%           D.vessel = true;
%           D.reserv = [false false true true true];
%           D.sol = true;
%           D.balance = 0.78; % balance the offset.
%           D.label = true;
% 
% output geometry struct with reshaped/postcalculated terms
% plot in current figure
% 
% take over axis of the plot if provided
%set('CurrentAxis',axis)
D.vessel = true;
D.reserv = [false false true true true];
D.sol = true;
D.balance = 0.78; % balance the offset.
D.label = true;
D.solgrid = false; % plot solgrid as well
P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end

% reshape normal vector
normal = [indiv.grid.nr_cc; indiv.grid.nz_cc];
geometry.normal_vector = normal;

balance = P.balance; %  0.78; 
[ir, iz] = offset_curve(indiv.grid.r_cc,indiv.grid.z_cc,-indiv.grid.sol_width_pol*(1-balance),normal);
[or, oz] = offset_curve(indiv.grid.r_cc,indiv.grid.z_cc,indiv.grid.sol_width_pol*balance,normal);

% reshape vessel
i_baffle = indiv.grid.i_baffle;
X_core_SOL = indiv.physics.x_core_sol;
mod(indiv.grid.vesrz,3);
Nvs = length(indiv.grid.vesrz)/3;
vesrz = reshape(indiv.grid.vesrz,Nvs,3);
[domain_bound,~] = find(vesrz(:,3)==1);
% construct arrays for neutral reservoirs
         clear romrz
         romrz{1}(1:2,:) = [vesrz(domain_bound(1):domain_bound(2),1:2); [ir(1:max(i_baffle(1),2))'  iz(1:max(i_baffle(1),2))']; vesrz(domain_bound(1),1:2)]';
         romrz{2}(1:2,:) = [vesrz(domain_bound(4):domain_bound(3),1:2); [or(1:max(i_baffle(1),2))'  oz(1:max(i_baffle(1),2))']; vesrz(domain_bound(4),1:2) ]';
         if  X_core_SOL  == 0  % only from stagnation point
         romrz{3}(1:2,:) = [flip(vesrz(domain_bound(5):domain_bound(6),1:2)); [or(i_baffle(1):i_baffle(2))' oz(i_baffle(1):i_baffle(2))']; vesrz(domain_bound(6),1:2)]';
         else % from baf to baf
         romrz{3}(1:2,:) = [flip(vesrz(domain_bound(4):domain_bound(6),1:2)); [or(i_baffle(1):i_baffle(2))' oz(i_baffle(1):i_baffle(2))']; vesrz(domain_bound(6),1:2)]';
         end
         romrz{4}(1:2,:) = [flip(vesrz(domain_bound(6):domain_bound(7),1:2)); [or(i_baffle(2):end)' oz(i_baffle(2):end)']; vesrz(domain_bound(7),1:2)]';
         romrz{5}(1:2,:) = [vesrz(domain_bound(8):domain_bound(9),1:2); [ir(i_baffle(2):end)' iz(i_baffle(2):end)']; vesrz(domain_bound(8),1:2)]';
geometry.romrz = romrz;
geometry.domain_bound = domain_bound;
geometry.vesrz = vesrz;
geometry.ir = ir;
geometry.or = or;
geometry.iz = iz;
geometry.oz = oz;

if P.vessel
plot(vesrz(:,1), vesrz(:,2),'k','LineWidth',2); hold on;
end

if P.sol 
    if P.solgrid
    plot([geometry.or(:)'; geometry.ir(:)' ],[geometry.oz(:)'; geometry.iz(:)'],'b','LineWidth',0.1); hold on;
    end
plot([geometry.or(:) geometry.ir(:) ],[geometry.oz(:) geometry.iz(:)],'b','LineWidth',1); hold on;
% plot(geometry.ir,geometry.iz,'b','LineWidth',1); hold on;
% plot(geometry.or,geometry.oz,'b','LineWidth',1); 
plot([geometry.or(1) geometry.ir(1) ],[geometry.oz(1) geometry.iz(1)],'b','LineWidth',1); 
plot([geometry.or(end) geometry.ir(end) ],[geometry.oz(end) geometry.iz(end)],'b','LineWidth',1); 
end

for irom = 1:5
    if P.reserv(irom)
        plot(romrz{irom}(1,:),romrz{irom}(2,:),'LineWidth',0.2); hold on;
    end
end

axis equal; xlabel('r [m]'); ylabel('z [m]'); hold off;
end

