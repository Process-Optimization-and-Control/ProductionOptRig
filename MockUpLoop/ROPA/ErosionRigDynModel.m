function [F,S_xx,S_zz,S_xz,S_xp,S_zp,x_var,z_var,u_var,p_var,diff,alg,L] = ErosionRigDynModel(par)
%    Creates a dynamic model of the rig and computes the sensitivity
%    matrices for EKF

% Inputs:
%    par = system parameters
%
% Outputs:
%   F: system integrator
%   S's: system sensitivities
%   x_var,z_var,u_var,p_var, diff,alg,L: Model in CasADi form 

% Other m-files required: none
% MAT-files required: none

% addpath('C:\Users\lab\Documents\casadi-windows-matlabR2016a-v3.4.5')                
import casadi.*

%% Parameters
%number of wells
n_w = par.n_w; %[]
%gas constant
R = par.R; %[m3 Pa K^-1 mol^-1]
%air molecular weigth
Mg = par.Mw; %[kg/mol] -- Attention: this unit is not usual

%properties
%density of water - dim:  nwells x 1
rho_l = par.rho_o; %[kg/m3]
%mixture viscosity
mu_mix = par.mu_oil;% [Pa s] 

%project - dim:  nwells x 1
% well length
L_w = par.L_w; %[m]
% well pipes cross section area
A_w = par.A_w;%[m2]

% riser length
L_r = par.L_r; %[m]
% riser pipes cross section area
A_r = par.A_r;%[m2]
%riser height 
H_r = par.H_r; %[m]
%well below injection
D = par.D_bh; %[m]

%% System states
% differential
%gas holdup 
m_g = MX.sym('m_g',n_w);            % 1:3 [1e-4 kg]
%water holdup 
m_l = MX.sym('m_l',n_w);            % 4:6 [kg]

% algebraic
%water rate from reservoir
w_l = MX.sym('w_l',n_w);            % 1:3 [1e-2 kg/s]
%total well production rate
w_total = MX.sym('w_total',n_w);    % 4:6 [1e-2 kg/s]

%riser head pressure
p_rh = MX.sym('p_rh',n_w);          % 7:9 [bar]
%pressure - before injection point (bottom hole)
p_bi = MX.sym('p_bi',n_w);          % 10:12 [bar]

%mixture density in system
rho_mix =  MX.sym('rho_mix',n_w);   % 13:15 [1e2 kg/m3]
%density gas
rho_g= MX.sym('rho_g',n_w);         % 16:18 [kg/m3]

%well outlet flowrate (gas)
w_gout = MX.sym('w_gout',n_w);      % 19:21 [1e-5 kg/s]
%riser head gas production rate gas
w_lout = MX.sym('w_lout',n_w);      % 22:24 [1e-2 kg/s]

%% System input
%gas lift rate
Q_gl = MX.sym('Q_gl',n_w);        % 1:3 [sL/min]
%valve oppening
vo = MX.sym('vo',n_w);            % 4:6 [0-1]
%pump outlet pressure
Ppump = MX.sym('Ppump',1);        % 7 [bar]

%% parameters
%%%%%%%%%
% fixed %
%%%%%%%%%
%room temperature
T = MX.sym('T',1); %[K]
%atmospheric pressure
p_atm = MX.sym('p_atm',1); %[bar]

%time transformation: CASADI integrates always from 0 to 1 and the USER does the time
%scaling with T --> sampling time
t_samp = MX.sym('t_samp',1); %[s]

% estimable
%scaled reservoir valve parameters
res_theta = MX.sym('res_theta',n_w);
%scaled top valve parameters
top_theta = MX.sym('top_theta',n_w);

%% Modeling
% Algebraic
%conversion
CR = 60*10^3; % [L/min] -> [m3/s] 
%reservoir outflow
f1 = -Ppump*ones(n_w,1)*1e5 + (w_l.*1e-2).^2.*(res_theta.*1e9)./(vo.^2.*rho_l) + p_bi.*1e5 ; 
% total system production 
f2 = - (w_total.*1e-2) + ((w_gout.*1e-5) + (w_lout.*1e-2));
%riser head pressure
f3 = -p_rh.*1e5 + (w_total.*1e-2).^2.*(top_theta.*1e8)./(rho_mix.*1e2) + p_atm.*1e5 ;
%before injection pressure 
f4 = -p_bi.*1e5 + (p_rh.*1e5 + (rho_mix.*1e2).*9.81.*H_r + 128.*mu_mix.*(L_w+L_r).*(w_l.*1e-2)./(3.14.*D.^4.*(rho_mix.*1e2)));
%mixture density
f5 = -(rho_mix.*1e2) + (((m_g.*1e-4) + m_l).* p_bi.*1e5.*Mg.*rho_l)./(m_l.*p_bi.*1e5.*Mg + rho_l.*R.*T.*(m_g.*1e-4));
%gas density (ideal gas law)
f6 = -rho_g + p_bi.*1e5.*Mg/(R*T);

% Simplifying assumption!
% liquid fraction in the mixture
xL = (m_l./((m_g.*1e-4) + m_l));
%Liquid outlet flowrate
f7 = -(w_lout.*1e-2) + xL.*(w_total.*1e-2);

% Total volume constraint
f8 = -(A_w.*L_w + A_r.*L_r) + (m_l./rho_l + (m_g.*1e-4)./rho_g);

% Differential
% gas mass balance
df1=  1e4*(-(w_gout.*1e-5) + Q_gl.*rho_g/CR);
% liquid mass balance
df2= -(w_lout.*1e-2) + (w_l.*1e-2); 

% Form the DAE system
diff = vertcat(df1,df2);
alg = vertcat(f1,f2,f3,f4,f5,f6,f7,f8);

% give fixed parameter values
alg = substitute(alg,p_atm,par.p_s);
alg = substitute(alg,T,par.T_r);

% concatenate the differential and algebraic states
x_var = vertcat(m_g,m_l);
z_var = vertcat(w_l,w_total,p_rh,p_bi,rho_mix,rho_g,w_gout,w_lout);
u_var = vertcat(Q_gl,vo,Ppump);
p_var = vertcat(res_theta,top_theta,t_samp);

%objective function 
L = 20*((w_lout(1)*1e-2)*CR/rho_l(1)) + 10*((w_lout(2)*1e-2)*CR/rho_l(2)) + 30*((w_lout(3)*1e-2)*CR/rho_l(3));
%end modeling

%% Casadi commands
%declaring function in standard DAE form (scaled time)
dae = struct('x',x_var,'z',z_var,'p',vertcat(u_var,p_var),'ode',t_samp*diff,'alg',alg);

%calling the integrator, the necessary inputs are: label; integrator; function with IO scheme of a DAE (formalized); struct (options)
F = integrator('F','idas',dae);

% ================================================
%       Calculating sensitivity matrices
% ================================================

S_xx = F.factory('sensStaStates',{'x0','z0','p'},{'jac:xf:x0'});
S_zz = F.factory('sensStaStates',{'x0','z0','p'},{'jac:zf:z0'});
S_xz = F.factory('sensStaStates',{'x0','z0','p'},{'jac:xf:z0'});

S_xp = F.factory('sensParStates',{'x0','z0','p'},{'jac:xf:p'});
S_zp = F.factory('sensParStates',{'x0','z0','p'},{'jac:zf:p'});

end
