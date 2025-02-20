function out = test_2PM_KR(f_pwr,f_mom,n_X,n_L,T_X,T_L,q_X,L,varargin)
% out = test_2PM_KR(f_pwr,f_mom,n_X,n_L,T_X,T_L,q_X,L,varargin)
% Use to evaluate Two-Point Model (2PM) preditions between X-point and target
% 
% For Kotov and Reiter 2PM (from Stangeby 2018)
%   - predict target temperature and density
% For Westermann and Westerhof 2PM (from DIV1D manual)
%   - predict target temperature and density
%   - predict upstream temperature
%
%
% INPUT:
%   - f_pwr  power loss fraction    
%   - f_mom  momentum loss fraction
%   - n_X    X-point density
%   - n_L    target density
%   - T_X    X-point temperature
%   - T_L    target temperature
%   - q_X    X-point heat flux
%   - L      Connection length from target to X-point
%   - Spu   % upstream pressure for K-R
% PARAMETERS
%   'print'         ==0 no print ==1 print comparison
%   'gamma'         ==6 default sheath heat transmission coeff.
%   'f_cnv'         ==0 default no upstream convective contribution.
%   'f_exp'         ==1 default no flux expansion considered = B_L / B_X  - 1
%   'pu'            ==-1 default, set to value if you want KR target density predictions 
% OUTPUT:
%   out = struct; {Wester,Kotov-Reiter,input}
%   out.T_u = [T_XW, T_XKR,   T_X];    % temperatures upstream 
%   out.T_t = [T_LW, T_LKR,   T_L];    % temperatures target
%   out.n_t = [n_LW, n_LKR,   n_L];    % densities target
%   out.f_cnv % X-point convective heat fraction
%
% NOTE the upstream temperature is predicted by Wester but an input for KR

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

%% Handling parameters
q_X = abs(q_X);
n_X = abs(n_X);
% Default paramters
D.print = '';
D.gamma = 6;
D.f_cnv = 0;
D.f_exp    = 1;
D.pu       = -1;
P = struct();
% Overwriting parameters
for k = 1:2:length(varargin), P.(varargin{k}) = varargin{k+1}; end
for k = fieldnames(D)'
 if ~isfield(P,k{1}), P.(k{1}) = D.(k{1}); end
end

% Adjust power and momentum loss fractions (cannot exceed 1)
f_pwr=max(min(0.9999,f_pwr),0);
f_mom=max(min(0.9999,f_mom),0);

%% Constants
e_charge = 1.602*1e-19;
m = 1.6605390401e-27;
kappa_0 = 2000;         % heat conductivity at 1 eV
mass    = 2*m;          
gamma_heat = P.gamma;   % sheath heat transmission factor

e_f = abs(P.f_exp-1);
% approximate flux expansion integral
f_exp_int = 1 - e_f/2 + e_f^2/3 - e_f^3/4 + e_f^4/5 - e_f^5/6;
if P.f_exp < 1
  %  disp('Te flux expansion might be interpreted as contraction') 
end
%% 2PM predictions equation
f_cnv_ = P.f_cnv;
count = 1;
diff = 1.0;
T_LW   = 0.0;
% line search with variable step size
while diff > 0.001 && count < 10000 % Armijo condition [nocedal 2008]
    T_XW    = (T_LW^(7/2) + 7.0 * q_X*(1-f_cnv_) * L * f_exp_int / 2.0 / kappa_0)^(2/7);
    T_LWnew = (mass / e_charge) * 2 * (q_X*(1.0-f_pwr))^2 / ( T_XW*(1-f_mom) )^2 / (gamma_heat * e_charge * n_X)^2 / (1 + e_f)^2;
    n_LW    = (n_X^3 / (q_X*(1.0-f_pwr))^2) * ( (T_XW*(1-f_mom))^3 * gamma_heat^2* e_charge^3 / 4 / mass ) * (1 + e_f)^2;
    diff      = abs(T_LWnew - T_LW);
    T_LW    = 0.1*T_LWnew+0.9*T_LW;
    if isnan(T_LW )
       disp('nan') 
    end
    count = count + 1;
    if count == 10000-2 
    disp('error: 2PM did not converge') 
    end
 end

%% Kotov and Reiter predictions
T_XKR    = T_X;
T_LKR = (mass / e_charge) * 2 * (q_X*(1.0-f_pwr))^2 / ( T_XKR*(1-f_mom) )^2 / (gamma_heat * e_charge * n_X)^2 / P.f_exp^2;
if P.pu == -1 
    n_LKR    =   gamma_heat^2 * (2*n_X*T_XKR*e_charge)^3 * (1-f_mom)^3 * P.f_exp^2 / 32 / mass / (q_X*(1.0-f_pwr))^2 ;
else % with SOLPS upstream pressure static + dynamic
    n_LKR    =   gamma_heat^2 * (P.pu)^3 * (1-f_mom)^3 * P.f_exp^2 / 32 / mass / (q_X*(1.0-f_pwr))^2 ;
end


%% print output
if length(P.print)>2
fprintf('%s \n','Comparison to Wester Formatted 2 Point Model')
fprintf('%s \n','inputs:')
fprintf('%s %7.1e %s %7.1e %s %4.1f %s \n','        upstream density n_X =', n_X, 'm^-3, upstream heat flux q_\parallel,X =', q_X, 'W/m^2, length of divertor leg L =', L, 'm')
fprintf('%s \n','results:')
fprintf('%s %8.2f %s %8.2f %s %8.2f %s \n','      momentum loss fraction f_mom =', f_mom,  ',   energy loss fraction f_pwr =', f_pwr, ' fraction convective', P.f_cnv,' ')
fprintf('%s %8.2f %s %s %s %8.2f %s \n','        Xpoint temperature 2PM T_X =', T_XW, 'eV,   Xpoint temperature',P.print ,'=', T_X, 'eV')
fprintf('%s %8.2f %s %s %s %8.2f %s \n','        target temperature 2PM T_L =', T_LW, 'eV,   target temperature',P.print ,'=', T_L, 'eV')
fprintf('%s %8.2e %s %s %s %8.2e %s \n','        target density 2PM     n_L =', n_LW, 'm^-3, target density    ',P.print ,'=', n_L, 'm^-3')
end
if length(P.print)>2
fprintf('%s \n','Comparison to Kotov-Reiter Formatted 2 Point Model')
fprintf('%s %8.2f %s %s %s %8.2f %s \n','        target temperature K-R T_L =', T_LKR, 'eV,   target temperature ',P.print ,'=', T_L, 'eV')
fprintf('%s %8.2e %s %s %s %8.2e %s \n','        target density K-R     n_L =', n_LKR, 'm^-3, target density ',P.print ,'=', n_L, 'm^-3')
end

out = struct; %{Wester,Kotov-Reiter,inputs}
out.T_u = [T_XW, T_XKR,   T_X];    % Temperatures upstream
out.T_t = [T_LW, T_LKR,   T_L];    % Temperatures target
out.n_t = [n_LW, n_LKR,   n_L];    % densities target
out.f_cnv = P.f_cnv; % convective X-point heat fraction used by Wester 2PM
out.pu = P.pu;
end