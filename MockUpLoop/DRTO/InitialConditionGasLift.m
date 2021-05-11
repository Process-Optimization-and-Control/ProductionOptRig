%Initial Condition Simplified
function [dx0,z0,u0,theta0] = InitialConditionGasLift(par)
%   Defines the initial state
%   pre-computed from experimental data using a SS estimation routine

% Inputs:
%    par = system parameters
%
% Outputs:
%   initial value of differential states (x)
%                    algebraic states (z)
%                    inputs (u)
%                    parameters (p)

% Other m-files required: none
% MAT-files required: InitialState_2021-03-03_104459_NoRTO_test2_1.mat

% previously computed
temp = load('InitialState_2021-03-03_104459_NoRTO_test2_1');

dx0 = vertcat(temp.m_g_0,temp.m_o_0);
z0 = vertcat(temp.w_ro_0,temp.w_pr_0,temp.p_rh_0,temp.p_bh_0,temp.rho_r_0,temp.rho_gr_0,temp.w_gr_0,temp.w_lr_0);
u0 = vertcat(temp.Q_gl_0,temp.vo_0,temp.vpump_0);
theta0  = vertcat(temp.res_theta_0,temp.val_theta_0);