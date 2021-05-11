function  par = ParametersGasLiftModel
%   Defines the system parameter that are fixed. The tuning parameters of
%   the three methods (SSRTO, ROPA and DRTO are also defined here

% Outputs:
%   par = system parameters

% Other m-files required: InitialConditionGasLift.m
% MAT-files required: none

%N.B.: the variable names do not match with the names in the paper. They
%match only in the model. However, the description matches
% Also, not all the parameters here are used in the simulations
%% 
%number of wells
par.n_w = 3; %[]
%gas constant
par.R = 8.314; %[m3 Pa/(K mol)]
%molecular weigth
par.Mw = 0.029; %[kg/mol] -- Attention: this unit is not usual

%% Properties
%density of water
par.rho_o = 996.57*ones(par.n_w,1); %[kg/m3]
%water viscosity
par.mu_oil = 1*0.000853*ones(par.n_w,1); %[Pa s or kg/(m s)]
%Gas densitiy at normal conditions of temperature and pressure
par.rho_g = 1.2041*ones(par.n_w,1); %[kg/m3]
%riser temperature
par.T_r = 23+273; %[K]
%atmospheric pressure
par.p_s = 1.01; %[bar]

%% Project
%well parameters - dim:  nwells x 1
%length
par.L_w = 1.8*ones(par.n_w,1); %[m]
%height
par.H_w  = 0*ones(par.n_w,1); %[m]
%diameter
par.D_w = 0.02*ones(par.n_w,1); %[m]
%well transversal area
par.A_w = pi.*(par.D_w/2).^2;%[m2]

%well before injection - [m]
par.L_bh = 0.4*ones(par.n_w,1);
par.H_bh = 0*ones(par.n_w,1);
par.D_bh = 0.02*ones(par.n_w,1);
par.A_bh = pi.*(par.D_bh/2).^2;%[m2]

%riser - [m]
par.L_r = 2.2*ones(par.n_w,1);
par.H_r = 2.2*ones(par.n_w,1);
par.D_r = 0.02*ones(par.n_w,1);
%riser areas
par.A_r = pi.*(par.D_r/2).^2;%[m2]

%% Upper limits 
%Max wellhead gas production rate 
par.QgMax = 7.5; % [sL/min] 

% Sampling time of the production optimization methods
par.T = 10; %[s]

%% Parameter Estimation
par.Sigma = eye(6); % since cov of P is super small, it is better to use the non-weigthed version

%% DRTO
[dx0,z0,u0,~] = InitialConditionGasLift(par);
% input dimension
par.nu = length(u0) - 4;
% differential state dimension
par.nx = length(dx0);
% algebraic state dimension
par.nz = length(z0);

% control horizon
% par.np = 6; % 60 s
par.np = 4; % 40 s
% prediction horizon
par.nm = par.np;

% input movement penalization weight
% par.RR = 10*eye(3,3);%diag([12 12 12]);
par.RR = 1*eye(3,3);%diag([12 12 12]);
% par.RR = 0*eye(3,3);%diag([12 12 12]);

% input movement constraint
par.dumax = 2.0; % [sL/min]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % testing the weight effect %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %conversion
% CR = 60*10^3; % [L/min] -> [m3/s]
% Ltest_term1 = 10*((z0(1)*1e-2)*CR/par.rho_o(1)) + 20*((z0(2)*1e-2)*CR/par.rho_o(2)) + 30*((z0(3)*1e-2)*CR/par.rho_o(3));
% Ltest_term2 = 0.5*(par.dumax*ones(3,1))'*par.RR*(par.dumax*ones(3,1));

%% for othogonal collocation
% Degree of interpolating polynomial
par.d = 3;

% Get collocation points
tau_root = [0 0.155051025721682 0.644948974278318 1];

% Coefficients of the collocation equation
par.C = zeros(par.d+1,par.d+1);

% Coefficients of the continuity equation
par.D = zeros(par.d+1, 1);

% Coefficients of the quadrature function
par.B = zeros(par.d+1, 1);

% Construct polynomial basis
for j=1:par.d+1
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:par.d+1
    if r ~= j
      coeff = conv(coeff, [1, -tau_root(r)]);
      coeff = coeff / (tau_root(j)-tau_root(r));
    end
  end
  % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  par.D(j) = polyval(coeff, 1.0);

  % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);
  for r=1:par.d+1
    par.C(j,r) = polyval(pder, tau_root(r));
  end

  % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);
  par.B(j) = polyval(pint, 1.0);
end



