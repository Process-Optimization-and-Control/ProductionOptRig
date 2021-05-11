% Main program
% Run Initialization file first

%%%%%%%%%%%%%%%%
% Get Variables
%%%%%%%%%%%%%%%%
% disturbances 
%valve opening [%] 
cv101 = P_vector(1);
cv102 = P_vector(2);
cv103 = P_vector(3);
% if value is [A] from 0.004 to 0.020
% if you want to convert to 0 (fully closed) to 1 (fully open)
% vo_n = (vo - 0.004)./(0.02 - 0.004);

%pump rotation [%]
pRate = P_vector(4);
% if value is [A] from 0.004 to 0.020
% if you want to convert to (min speed - max speed)
% goes from 12% of the max speed to 92% of the max speed
% pRate = 12 + (92 - 12)*(P_vector(4) - 0.004)./(0.02 - 0.004);

% always maintain the inputs greater than 0.5
% inputs computed in the previous MPC iteration
% Note that the inputs are the setpoints to the gas flowrate PID's
fic104sp = P_vector(5);
fic105sp = P_vector(6); 
fic106sp = P_vector(7);
%current inputs of the plant
u0old=[P_vector(5);P_vector(6);P_vector(7)];


% cropping the data vector
nd = size(I_vector,2);
dataCrop = (nd - BufferLength + 1):nd;

% liquid flowrates [L/min]
fi101 = I_vector(1,dataCrop);
fi102 = I_vector(2,dataCrop);
fi103 = I_vector(3,dataCrop);

% actual gas flowrates [sL/min]
fic104 = I_vector(4,dataCrop);
fic105 = I_vector(5,dataCrop);
fic106 = I_vector(6,dataCrop);

% pressure @ injection point [mbar g]
pi105 = I_vector(7,dataCrop);
pi106 = I_vector(8,dataCrop);
pi107 = I_vector(9,dataCrop);

% reservoir outlet temperature [oC]
ti101 = I_vector(10,dataCrop);
ti102 = I_vector(11,dataCrop);
ti103 = I_vector(12,dataCrop);

% DP @ erosion boxes [mbar]
dp101 = I_vector(13,dataCrop);
dp102 = I_vector(14,dataCrop);
dp103 = I_vector(15,dataCrop);

% top pressure [mbar g]
% for conversion [bar a]-->[mbar g]
% ptop_n = ptop*10^-3 + 1.01325;
pi101 = I_vector(16,dataCrop);
pi102 = I_vector(17,dataCrop);
pi103 = I_vector(18,dataCrop);

% reservoir pressure [bar g]
% for conversion [bar g]-->[bar a]
% ptop_n = ptop + 1.01325;
pi104 = I_vector(19,dataCrop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE GOES HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of measurements in the data window
dss = size(pi104,2);

% we do no consider SS detection here
flagSS = 1;

% check for first iteration
if ~exist('dxHatkk','var') 
    %initial condition
    [dx0,z0,u0,theta0] = InitialConditionGasLift(par);

    dxHatkk = dx0;
    zHatkk = z0;
    thetakk = theta0;

    load('EKFconf');   
    qThres = ekf.qThresk;
    Pkk = ekf.Pkk;
    ekf.R = noise.output; %added (different from Rig Implementation)
end

if flagSS == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimating SS model parameters (dynamic) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yPlant = [fi101; 
              fi102;
              fi103;
              1.01325 + 10^-3*pi101; %[mbarg]-->[bar a];
              1.01325 + 10^-3*pi102;
              1.01325 + 10^-3*pi103];
          
    uPlant = [fic104; %conversion [L/min] --> [kg/s]
              fic105;
              fic106;
              cv101*ones(1,dss); %workaround - i just have the last measurement here. Since it is the disturbance, it doesn't really matter;
              cv102*ones(1,dss);
              cv103*ones(1,dss);
              pi104 + 1.01325];            
    
          try
              % running casadi - getting the last measurement and last
              % input that generated that measuremnt
              [dxHatkk,zHatkk,thetakk,Pkk,qThres] = RExtendedKalmanFilter(yPlant(:,end),uPlant(:,end - 10),dxHatkk,zHatkk,thetakk,Pkk,F_model,S_xx,S_zz,S_xz,S_xp,S_zp,ekf,par);
              
              % everything normal
              flagEst = 1;
          catch
              warning('Dynamic Estimation Problem!');
              beep
              % estimation problem
              flagEst = 0;
              
              % Note that in this case, we dont update
              % dxHatkk,zHatkk,thetakk,Pkk,qThres
          end
    
    
    if flagEst == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optimizing SS model parameters %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getting mean value of the inputs as initial guess for optimization
        uGuess = mean(uPlant,2);
        %uGuess(1:3) = [2.5; 2.5; 2.5];
        
        [uOpt,xOptHat,zOptHat,phi,flagOpt] = ErosionRigSSOptimization(dxHatkk,zHatkk,thetakk,uGuess,par);

        if flagOpt == 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Send Variable
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % compute new values for the gas flow rate setpoints
            %Filter for optimum
            fic104sp = u0old(1) + OptConf.ku.*(uOpt(1) - u0old(1));
            fic105sp = u0old(2) + OptConf.ku.*(uOpt(2) - u0old(2));
            fic106sp = u0old(3) + OptConf.ku.*(uOpt(3) - u0old(3));

            O_vector = vertcat(fic104sp,fic105sp,fic106sp)';
            SS = 1;
            Estimation = 1;
            Optimization = 1;
            Result = phi; 
            Parameter_Estimation = thetakk';
            State_Variables_Estimation = (par.H*zHatkk)';
            State_Variables_Optimization = (par.H*zOptHat)';
            Optimized_Air_Injection = uOpt';
        else
            %%%%%%%%%
            %(dummy)%
            %%%%%%%%%
            % compute new values for the gas flow rate setpoints
            O_vector = vertcat(fic104sp,fic105sp,fic106sp)';

            SS = 1;
            Estimation = 1;
            Optimization = 0;
            Result = 0;
            Parameter_Estimation = [0,0,0,0,0,0];
            State_Variables_Estimation = [0,0,0,0,0,0];
            State_Variables_Optimization = [0,0,0,0,0,0];
            Optimized_Air_Injection = [0,0,0];
        end
        
    else  
        %%%%%%%%%
        %(dummy)%
        %%%%%%%%%
        % compute new values for the gas flow rate setpoints
        O_vector = vertcat(fic104sp,fic105sp,fic106sp)';

        SS = 1;
        Estimation = 0;
        Optimization = 0;
        Result = 0;
        Parameter_Estimation = [0,0,0,0,0,0];
        State_Variables_Estimation = [0,0,0,0,0,0];
        State_Variables_Optimization = [0,0,0,0,0,0];
        Optimized_Air_Injection = [0,0,0];
    end
    
    
else
    %%%%%%%%%
    %(dummy)%
    %%%%%%%%%
    % compute new values for the gas flow rate setpoints
    O_vector = vertcat(fic104sp,fic105sp,fic106sp)';
    
    SS = 0;
    Estimation = 0;
    Optimization = 0;
    Result = 0;
    Parameter_Estimation = [0,0,0,0,0,0];
    State_Variables_Estimation = [0,0,0,0,0,0];
    State_Variables_Optimization = [0,0,0,0,0,0];
    Optimized_Air_Injection = [0,0,0];
    
end
