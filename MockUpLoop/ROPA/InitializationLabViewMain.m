%clear
%clc

%% Model tuning
%paramters
par = ParametersGasLiftModel;
%initial condition
[dx0,z0,u0,theta0] = InitialConditionGasLift(par);
% model integrators
[F_model,S_xx,S_zz,S_xz,S_xp,S_zp,x_var,z_var,u_var,p_var,difEq,algEq,L] = ErosionRigDynModel(par);

%states to measurement mapping function
par.nMeas = 6;
H = zeros(6,24);
H(1,1) = 1e-2*60*1e3/par.rho_o(1); %wro-oil rate from reservoir, well 1 [1e2 kg/s] --> [1e2 L/min]
H(2,2) = 1e-2*60*1e3/par.rho_o(2); %wro-oil rate from reservoir, well 2
H(3,3) = 1e-2*60*1e3/par.rho_o(3); %wro-oil rate from reservoir, well 3
H(4,7) = 1; %prh - riser head pressure well 1
H(5,8) = 1; %prh - riser head pressure well 2
H(6,9) = 1; %prh - riser head pressure well 3
par.H = H;
par.Hext = [zeros(6,6), par.H, zeros(6,6)]; %differential states, algebraic states, parameters

%% Optimization tuning
OptConf.ku = 0.4; % input filter

%% Initializing filter
% initializing parameter estimate array
dxk = dx0;
zk = z0;

% %Filter - for reduced EKF. Not used here
qThresk = 36; 