function [uOpt,xHat,zHat,phi,flagOpt] = ErosionRigDynOptimization(x0_hat,z0_hat,u_1,p_1,x_var,z_var,u_var,p_var,diff,alg,L,par)
%    Economic optimization of the dynamic rig model 

% Inputs:
%    x0_hat, z0_hat, p_1 = previously estimated states and parameters
%    u_1 = current input value
%    x_var,z_var,u_var,p_var,diff,alg,L = dynamic model equations in CasADi
%    par = system parameters
%
% Outputs:
%   uOpt = computed inputs
%   xHat, zHat = states computed at new SS optimum
%   phi = economic OF value at the optimum
%   flagOpt  = optimization didn't converge == 0 | success == 1

% Other m-files required: OptimizationBoundsGasLiftRiser2.m
% MAT-files required: none

    %addpath ('C:\Users\lab\Documents\casadi-windows-matlabR2016a-v3.4.5')
    import casadi.*
    
   %% Parameters  
    %Variable bounds
    [lbx,lbz,lbu,ubx,ubz,ubu,~,~] = OptimizationBoundsGasLiftRiser2(par);

    %% Model Equations - modifying objective function: regularization term
    % previous input value
    Uk_1 = MX.sym('Uk_1',par.nu);
    % current input value
    Uk = MX.sym('Uk',par.nu);

    % objective function for the controller (economic + input movement
    % constraint)
    Q =  -L + 1/2 * (Uk - Uk_1)'*par.RR*(Uk - Uk_1);
    
    % Function for the system model
    f = Function('f',{x_var,z_var,u_var,p_var,Uk,Uk_1},{diff,alg,Q});

   %% Controller
   % defining variables
    w = {};
    w0 = [];
    lbw =[];
    ubw = [];
    
    % constraints
    lbg = [];
    ubg = [];
    g = {};
    
    % initialize objective function
    J = 0;
    
    % "lifting" initial condition
    Xk = MX.sym('X0',par.nx);
    w = {w{:},Xk};
    lbw = [lbw,x0_hat];
    ubw = [ubw,x0_hat];
    w0 = [w0;x0_hat];
    
    % getting the initial input value
    Uprev = u_1(1:par.nu); 

%% Looping until timeend
% Formulate the NLP
for k=0:par.np-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],par.nu);
    w = {w{:}, Uk};
    lbw = [lbw; lbu(1:par.nu)];
    ubw = [ubw; ubu(1:par.nu)];
    w0 = [w0;  u_1(1:par.nu)];
    
    % constraint on max gas availability
    g = {g{:},sum(Uk)};
    lbg = [lbg;0];
    ubg = [ubg;par.QgMax];
    
    % creating input movement variables
    duk = Uk - Uprev;
    g = {g{:},duk};
    
    % checking the bounds (for control horizon ~= prediction horizon)
    if k > par.nm
        lbg = [lbg;zeros(par.nu,1)];
        ubg = [ubg;zeros(par.nu,1)];
    else
        lbg = [lbg;-par.dumax*ones(par.nu,1)];
        ubg = [ubg;par.dumax*ones(par.nu,1)];
    end

    % Creating states at collocation points
    Xkj = {};
    for j=1:par.d
        %differential states
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], par.nx);
        w = {w{:}, Xkj{j}};
        lbw = [lbw; lbx];
        ubw = [ubw; ubx];
        w0 = [w0; x0_hat];
        
        %algebraic states
        Zkj{j} = MX.sym(['Z_' num2str(k) '_' num2str(j)], par.nz);
        w = {w{:}, Zkj{j}};
        lbw = [lbw; lbz];
        ubw = [ubw; ubz];
        w0 = [w0; z0_hat];
    end

    % Loop over collocation points
    Xk_end = par.D(1)*Xk;
    for j=1:par.d
       % Expression for the state derivative at the collocation point
       xp = par.C(1,j+1)*Xk;
       for r=1:par.d
           xp = xp + par.C(r+1,j+1)*Xkj{r};
       end

       % Append collocation equations
        [fj,gj,qj] = f(Xkj{j},Zkj{j},vertcat(Uk,u_1(par.nu + 1:end)),p_1,Uk,Uprev);
       
       % dynamic states trajectories
       g = {g{:}, par.T*fj - xp};
       lbg = [lbg; zeros(par.nx,1)];
       ubg = [ubg; zeros(par.nx,1)];
       
       % enforcing algebraic constraints at every collocation point
       g = {g{:}, gj};
       lbg = [lbg; zeros(par.nz,1)];
       ubg = [ubg; zeros(par.nz,1)];
       
       % Add contribution to the end state
       Xk_end = Xk_end + par.D(j+1)*Xkj{j};

       % Add contribution to quadrature function
       J = J + par.B(j+1)*qj*par.T;
    end

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], par.nx);
    w = {w{:}, Xk};
    lbw = [lbw; lbx];
    ubw = [ubw; ubx];
    w0 = [w0; x0_hat];
    
    % Add equality constraint
    g = {g{:}, Xk_end-Xk};
    lbg = [lbg; zeros(par.nx,1)];
    ubg = [ubg; zeros(par.nx,1)];
    
    % updating Uprev
    Uprev = Uk;
end


%% Solving optimization problem
% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);

%% Extracting solution
w_opt = full(sol.x);
%c_opt = full(sol.g);
%lag_opt = full(sol.lam_g);

uOpt = w_opt(par.nx + 1:par.nx + par.nu);
xHat = w_opt(par.nx + par.nu + par.d*(par.nx + par.nz) + 1:par.nx + par.nu + par.d*(par.nx + par.nz) + par.nx); % next instant
zHat = w_opt(par.nx + par.nu + par.d*(par.nx + par.nz) + par.nx + 1:par.nx + par.nu + par.d*(par.nx + par.nz) + par.nx + par.nz);
phi = -full(sol.f);

if solver.stats.success == 1
    flagOpt = 1;
else
    flagOpt = 0;
    msg = ['Error optimizing model'];
    error(msg);
end


end

