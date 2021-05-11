function [lbx,lbz,lbu,ubx,ubz,ubu,lbtheta,ubtheta] = OptimizationBoundsGasLiftRiser2(par)
%   Defines the lower and upper bounds for the model variables

% Inputs:
%    par = system parameters
%
% Outputs:
%   bounds on differential states (x)
%             algebraic states (z)
%             inputs (u)
%             parameters (p)

% Other m-files required: none
% MAT-files required: none

%number of wells
n_w = par.n_w;

%% States
%gas holdup [kg]
m_g_lb = 0*ones(n_w,1);
m_g_ub = 1e1*ones(n_w,1);

%oil holdup [kg]
m_l_lb = 0.*ones(n_w,1);
m_l_ub = 1e3.*ones(n_w,1);

%water rate from reservoir [1e-2 kg/s]
w_l_lb = 0.*ones(n_w,1);
w_l_ub = 20.*ones(n_w,1);

% total production rate [1e-2 kg/s]
w_total_lb = 0.*ones(n_w,1);
w_total_ub = 60.*ones(n_w,1); 

%riser head pressure [bar]
p_rh_lb = 0.*ones(n_w,1);
p_rh_ub = 1.500.*ones(n_w,1); %49

%pressure - before injection point (bottom hole) [bar]
p_bi_lb = 0.*ones(n_w,1);
p_bi_ub = 2.*ones(n_w,1);

%mixture density in riser [1e2 kg/m3]
rho_mix_lb = 0.*ones(n_w,1);
rho_mix_ub = 15.*ones(n_w,1);

%density gas
rho_g_lb= 0.*ones(n_w,1);
rho_g_ub= 10.*ones(n_w,1);

% total production rate [1e-5 kg/s]
w_gout_lb = 0.*ones(n_w,1);
w_gout_ub = 15.*ones(n_w,1); 

% total production rate [1e-2 kg/s]
w_lout_lb = 0.*ones(n_w,1);
w_lout_ub = 60.*ones(n_w,1); 

%% Inputs
%gas lift rate [sL/min]
Q_gl_lb = 1.*ones(n_w,1); %0.5
Q_gl_ub = 5.5*ones(n_w,1); 

% valve opening [0-1]
vo_lb = 0.1*ones(n_w,1);
vo_ub = 0.995*ones(n_w,1); 

% pump rotation [bar]
Ppump_lb = 1.01325*ones(1,1); %atmosferic pressure
Ppump_ub = 2*1.01325.*ones(1,1); 

%% Parameters
%Theta
res_theta_lb= 0.*(ones(n_w,1));
res_theta_ub= 10*(ones(n_w,1));

%Riser valve caracteristics
top_theta_lb=0.*ones(n_w,1);
top_theta_ub=10.*ones(n_w,1);

%% 
lbx = vertcat(m_g_lb,m_l_lb);
lbz = vertcat(w_l_lb,w_total_lb,p_rh_lb,p_bi_lb,rho_mix_lb,rho_g_lb,w_gout_lb,w_lout_lb);
lbu = vertcat(Q_gl_lb,vo_lb,Ppump_lb);
lbtheta = vertcat(res_theta_lb,top_theta_lb);

ubx = vertcat(m_g_ub,m_l_ub);
ubz = vertcat(w_l_ub,w_total_ub,p_rh_ub,p_bi_ub,rho_mix_ub,rho_g_ub,w_gout_ub,w_lout_ub);
ubu = vertcat(Q_gl_ub,vo_ub,Ppump_ub);
ubtheta = vertcat(res_theta_ub,top_theta_ub);