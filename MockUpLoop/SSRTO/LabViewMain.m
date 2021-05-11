% Main program
% Run Initialization file first

%%%%%%%%%%%%%%%%
% Get Variables
%%%%%%%%%%%%%%%%
% disturbances 
%valve opening [-] 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOST RECENT VALUE IS THE LAST ONE! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% number of measurements in the data window
dss = size(pi104,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE GOES HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run SS identification
%relevant measurements are only the liquid flowrates
yPlantArray = [fi101;fi102;fi103];

[flagSS,pSS] = SSDetection(yPlantArray);

if flagSS == 1 % we are at steady state at the current instant
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimating SS model parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    % if the last measurement is considered at SS, we assume that the whole period between the current and past RTO execution is at SS (mean filtering)      
    uEst = mean(uPlant(:,end - 10:end),2);
    yEst = mean(yPlant(:,end - 10:end),2);

    % running casadi
    [thetaHat,xEstHat,zEstHat,yEstHat,~,flagEst] = ErosionRigSSEstimation(xEstHat,zEstHat,thetaHat,uEst,yEst,par);

    if flagEst == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optimizing SS model parameters %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uGuess = uEst;
        
        [uOpt,xOptHat,zOptHat,phi,flagOpt] = ErosionRigSSOptimization(xEstHat,zEstHat,thetaHat,uGuess,par);
        
        if flagOpt == 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Send Variable
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Filter for optimum
            fic104sp = u0old(1) + OptConf.ku.*(uOpt(1) - u0old(1));
            fic105sp = u0old(2) + OptConf.ku.*(uOpt(2) - u0old(2));
            fic106sp = u0old(3) + OptConf.ku.*(uOpt(3) - u0old(3));

            O_vector = vertcat(fic104sp,fic105sp,fic106sp)';
            SS = 1;
            Estimation = 1;
            Optimization = 1;
            Result = phi;
            Parameter_Estimation = thetaHat';
            State_Variables_Estimation = (par.H*zEstHat)';
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
