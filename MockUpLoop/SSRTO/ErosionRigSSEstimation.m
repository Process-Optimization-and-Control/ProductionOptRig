function [thetaHat,xHat,zHat,yHat,Psi,flagEst] = ErosionRigSSEstimation(xGuess,zGuess,thetaGuess,uk,yk,par)
%    Estimates Parameters of the steady-state rig model 

% Inputs:
%    xGuess, zGuess = guess for the states 
%    thetaGuess = parameters guess
%    uk = current inputs (measured)
%    yk = current measurements
%    par = system parameters
%
% Outputs:
%   thetaHat = estimated parameters
%   xHat, zHat, yHat = estimated states and model prediction
%   Psi = estimation OF value
%   flagEst  = estimation didn't converge == 0 | success == 1

% Other m-files required: none
% MAT-files required: none

%addpath ('C:\Users\lab\Documents\casadi-windows-matlabR2016a-v3.4.5')                
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
D = par.D_bh;%[m]

%% System states
%N.B.: In the SS model, there is no algebraic and
%differential states. However, we kept differentiation the x's and z's for
%matching the dynamic model states. 

% differential
%gas holdup 
m_g = MX.sym('m_gr',n_w);               % 1:3 [1e-4 kg]
%oil holdup 
m_l = MX.sym('m_l',n_w);                % 4:6 [kg]

% algebraic
%oil rate from reservoir
w_l = MX.sym('w_l',n_w);                % 1:3 [1e-2 kg/s]
%riser head total production rate
w_total = MX.sym('w_total',n_w);        % 4:6 [1e-2 kg/s]

%riser head pressure
p_rh = MX.sym('p_rh',n_w);              % 7:9 [bar]
%pressure - below injection point (bottom hole)
p_bi = MX.sym('p_bi',n_w);              % 10:12 [bar]

%mixture density in riser
rho_mix =  MX.sym('rho_mix',n_w);       % 13:15 [100 kg/m3]
%density gas
rho_g= MX.sym('rho_g',n_w);             % 16:18 [kg/m3]

%riser head gas production rate gas
w_gout = MX.sym('w_gout',n_w);          % 19:21 [1e-5 kg/s]
%riser head gas production rate gas
w_lout = MX.sym('w_lout',n_w);          % 22:24 [1e-2 kg/s]

%% System input
%gas lift rate
Q_gl = MX.sym('Q_gl',n_w);      % 1:3 [L/min]
%valve oppening
vo = MX.sym('vo',n_w);          % 4:6 [0-1]
%pump outlet pressure
Ppump = MX.sym('Ppump',1);      % 7 [bar]

%% Parameters
% fixed
%room temperature
T = MX.sym('T',1); %[K]
%atmospheric pressure
p_atm = MX.sym('p_atm',1); %[bar]

% estimable
%scaled reservoir valve parameters
res_theta = MX.sym('res_theta',n_w);
%scaled top valve parameters
top_theta = MX.sym('top_theta',n_w);

%% Modeling
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
f5 = -(rho_mix.*1e2) + (((m_g.*1e-4) + m_l).*p_bi.*1e5.*Mg.*rho_l)./(m_l.*p_bi.*1e5.*Mg + rho_l.*R.*T.*(m_g.*1e-4));
%gas density (ideal gas law)
f6 = -rho_g + p_bi.*1e5.*Mg/(R*T);

% Simplifying assumption! -> used in a different form
% Gas mass holdup
f7 = -(m_g.*1e-4) + ((m_g.*1e-4) + m_l).* ((Q_gl./(CR./rho_g))./(Q_gl./(CR./rho_g) + (w_l.*1e-2))); 
% Oil mass holdup
f8 = -m_l + (rho_l).*((A_w.*L_w + A_r.*L_r) - (((m_g.*1e-4))./rho_g));

% Differential
% gas mass balance -> dropped the scaling factor
df1= -(w_gout.*1e-5) + Q_gl./(CR./rho_g);
% liquid mass balance
df2= -(w_lout.*1e-2) + (w_l.*1e-2);  

% Form system of equations
diff = vertcat(df1,df2);
alg = vertcat(f1,f2,f3,f4,f5,f6,f7,f8);

% give parameter values
alg = substitute(alg,p_atm,par.p_s);
alg = substitute(alg,T,par.T_r);

% concatenate the variables
x_var = vertcat(m_g,m_l);
z_var = vertcat(w_l,w_total,p_rh,p_bi,rho_mix,rho_g,w_gout,w_lout);
p_var = vertcat(Q_gl,vo,Ppump,res_theta,top_theta);

%end modeling

%% Estimation problem
% ===================================
%     Declaring 
% ===================================
% decision variables
w = {x_var,z_var,p_var};
% constraints (model)
g = {diff,alg};

%objective function
% computing residual
eps = (par.H*z_var - yk);

J = eps'*par.Sigma*eps;
%J = eps'*eps;

% formalize it into an NLP problem
nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

% Assign solver
options = struct;
options.ipopt.print_level = 5;
% options.ipopt.acceptable_constr_viol_tol = 0.01;
Est = nlpsol('solver','ipopt',nlp,options);

% ===================================
%     Solving
% ===================================
%Variable bounds
[lbx,lbz,~,ubx,ubz,~,lbtheta,ubtheta] = OptimizationBoundsGasLiftRiser2(par);

wk = [xGuess;zGuess;uk;thetaGuess];
lbw = [lbx;lbz;uk;lbtheta]; %lbtheta
ubw = [ubx;ubz;uk;ubtheta]; %ubtheta
% fixing uk

lbg = zeros(length(xGuess) + length(zGuess),1); %SS: diff and alg == 0
ubg = zeros(length(xGuess) + length(zGuess),1);

% Solve
sol = Est('x0',wk,'lbg',lbg,'ubg',ubg,'lbx',lbw,'ubx',ubw);

flagEst = Est.stats.success;
if Est.stats.success ~=1
    msg = ['Error estimating parameters, model'];
    error(msg);
end

%clc %clean ipopt output

% Extract Solution
xHat = full(sol.x(1:length(xGuess)));
zHat = full(sol.x(length(xGuess) + 1:length(xGuess) + length(zGuess)));
yHat = par.H*zHat;%conversion [kg/s] --> [L/min]
thetaHat = full(sol.x(length(xGuess) + length(zGuess) + length(uk) + 1:end));
Psi = full(sol.f);

end
