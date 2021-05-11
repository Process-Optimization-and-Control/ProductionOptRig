%clear
%clc

%% SS detection Configuration
SSConf.sspn = 5; %Number of points in SS to assure that system is at SS

%SS detection for the whole plant (MPA)
SSConf.dss = 10; %length of the moving window used for SS identification
SSConf.ny = 3;
SSConf.rThres = [1.5;2;2];
SSConf.lambda1 = [0.05;0.05;0.05];
SSConf.lambda2 = [0.1;0.1;0.1];
SSConf.lambda3 = [0.1;0.1;0.1];

%% Model tuning
%paramters
par = ParametersGasLiftModel;
%initial condition
[x0,z0,u0,theta0] = InitialConditionGasLift(par);

%states to measurement mapping function
par.nMeas = 6;
H = zeros(par.nMeas,length(z0));
H(1,1) = 1e-2*60*1e3/par.rho_o(1); %wro-oil rate from reservoir, well 1 [1e2 kg/s] --> [L/min]
H(2,2) = 1e-2*60*1e3/par.rho_o(2); %wro-oil rate from reservoir, well 2
H(3,3) = 1e-2*60*1e3/par.rho_o(3); %wro-oil rate from reservoir, well 3
H(4,7) = 1; %prh - riser head pressure well 1
H(5,8) = 1; %prh - riser head pressure well 2
H(6,9) = 1; %prh - riser head pressure well 2
par.H = H;

%% Optimization tuning
OptConf.ku = 0.4; % input filter

%% initializing estimation
xEstHat = x0;
zEstHat = z0;
thetaHat = theta0;

