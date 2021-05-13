% Runs a mock-up RTO problem, where the model and the system are equal and the disturbance setup is pre-determined. 
% The system setup is based on the paper:
%   Implementation of Steady-state Real-time Optimization Using Transient Measurements on Experimental Rig
%   by Matias et al.

% Other m-files required: SmoothPlantParameters2.mat; Disturbances2.mat; RigNoise.mat; InitialState_2021-03-03_104459_NoRTO_test2_1.mat
% MAT-files required: 
%       For SSRTO:
%           1. SSDetection.m
%           2. ErosionRigSSEstimation.m
%           3. ErosionRigSSOptimization.m

%       For Mockup Labview Interface:
%           1. InitializationLabViewMain.m
%           2. LabViewMain.m

%       For "Plant":
%           1. ErosionRigDynModel.m

%       Bounds, Parameters and Initial conditon:
%           1. InitialConditionGasLift.m
%           2. OptimizationBoundsGasLiftRiser2.m
%           3. ParametersGasLiftModel.m
clear
close all
clc

%saving data
name = 'SSRTO_DT_Cycle';
%noise seed
rng('default')

%% Loading .mat files
% previously defined disturbance profiles
disturbances = load('Disturbances3'); 

% previously computed system parameter profiles - used as "plant" true
% parameters
parProfile = load('SmoothPlantParameters3'); 

% rig noise characteristics - compute previously from actual rig data
noise = load('RigNoise'); 

%% Simulation tuning
% Buffer length
BufferLength = 60; %[s]

% Optimization sampling time
nExec = 10; %[s]

%plant parameters
parPlant = ParametersGasLiftModel;

%simulation parameters
nInit = 0; %[s]
nFinal = size(parProfile.thetaPlant,2); %[s] ! Adding a 60 second buffer in comparison to experimental duration (1244 ~ 22min) to emulate the measurement buffer of labview, which is "filled" before the experiment
parPlant.T = 1; %rig measurements samplint time[s]
tgrid = (nInit:parPlant.T:nFinal)/60; %[min] one measurements per second

%initial condition
[dxPlant0,zPlant0,uPlant0,thetaPlant0] = InitialConditionGasLift(parPlant);

%states to measurement mapping function
parPlant.nMeas = 6;
parPlant.H = zeros(6,24);
parPlant.H(1,1) = 1e-2*60*1e3/parPlant.rho_o(1); %wro-oil rate from reservoir, well 1 [1e-2 kg/s] --> [L/min]
parPlant.H(2,2) = 1e-2*60*1e3/parPlant.rho_o(2); %wro-oil rate from reservoir, well 2
parPlant.H(3,3) = 1e-2*60*1e3/parPlant.rho_o(3); %wro-oil rate from reservoir, well 3
parPlant.H(4,7) = 1; %prh - riser head pressure well 1
parPlant.H(5,8) = 1; %prh - riser head pressure well 2
parPlant.H(6,9) = 1; %prh - riser head pressure well 3

% Model representing the rig
[F,~,~,~,~,~,~,~,~,~,~,~,~] = ErosionRigDynModel(parPlant);

% mockup controller dynamic parameters 
tau_C = 1;

%% Run configuration file
InitializationLabViewMain %here we use the same syntax as in the rig

%% Initializing simulation
% Plant states
dxk = dxPlant0;
zk = zPlant0;

% Inputs
uk = [uPlant0(1);                   %FIC-104 [sl/min]
      uPlant0(2);                   %FIC-105 [sl/min]
      uPlant0(3);                   %FIC-106 [sl/min]
      disturbances.values(1,1);     %CV-101 opening [-]
      disturbances.values(2,1);     %CV-102 opening [-]
      disturbances.values(3,1);     %CV-103 opening [-]
      disturbances.values(4,1)];    %PI-104 [bar]

% Setpoints for the gas lift controllers
O_vector = [uk(1); uk(2); uk(3)];
  
% "Plant" parameters
thetak = thetaPlant0;

%% Run mock-up loop
% arrays for plotting
%%%%%%%%%%%%%%%%
% "Plant" Data %
%%%%%%%%%%%%%%%%
xPlantArray = dxk;
zPlantArray = zk;
uPlantArray = uk;
uSPPlantArray = O_vector; % decision variable of the optimization problem are the SP of the PI controllers that regulate the gas flowrate
measPlantArray = parPlant.H*zk;
thetaPlantArray = thetak;
ofPlantArray = 20*(measPlantArray(1)) + 10*(measPlantArray(2)) + 30*(measPlantArray(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Production Optimization Methods Data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we save the same data that is saved in the Rig's labview code
flagArray = []; % flag SSD, Model Adaptation and Economic Optimization [flag == 0 (failed), == 1 (success)]
ofArray = [];   % OF value computed bu the Economic Optimization
thetaHatArray = []; % estimated parameters
yEstArray = []; % model prediction with estimated states
yOptArray = []; % model prediction @ new optimum
uOptArray = []; % computed inputs (u_k^\star)
uImpArray = []; % filtered inputs to be implemented (u_{k+1})

for kk = 1:nFinal
    
    % printing the loop evolution in minutes
    fprintf('     kk >>> %6.4f [min]\n',tgrid(kk + 1))   
    
    % integrating the "plant" (the model sampling time is defined by parPlant.T)
    Fend = F('x0',dxk,'z0',zk,'p',[uk;thetak;parPlant.T]); 
    
    %extracting the results (from Casadi symbolic to numerical)
    dxk = full(Fend.xf);
    zk = full(Fend.zf);

    % saving the results
    xPlantArray = [xPlantArray, dxk];
    zPlantArray = [zPlantArray, zk];
    measPlantArray = [measPlantArray, parPlant.H*zk + noise.output*randn(6,1)]; %adding artificial noise to the measurements
    ofPlantArray = [ofPlantArray, 20*(measPlantArray(1,end)) + 10*(measPlantArray(2,end)) + 30*(measPlantArray(3,end));];
    
    % we execute the production optimization:
    % a. after the initial [BufferLength]-second buffer
    % b. every [nExec] seconds 
    if kk > BufferLength && rem(kk,nExec) == 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rearranging the data vectors and units here, so they can the be exactly %
        % the same as in the actual rig. The goal is that LabViewRTO.m can be     % 
        % directly plug in the Labview interface and it will work                 %      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measurement buffer (dim = 19 X BufferLength)
        I_vector = [measPlantArray(1,(kk - BufferLength + 2):kk + 1); % FI-101 [L/min]
                    measPlantArray(2,(kk - BufferLength + 2):kk + 1); % FI-102 [L/min]
                    measPlantArray(3,(kk - BufferLength + 2):kk + 1); % FI-103 [L/min]
                    uPlantArray(1,(kk - BufferLength + 1):kk);    % FI-104 [sL/min]
                    uPlantArray(2,(kk - BufferLength + 1):kk);    % FI-105 [sL/min]
                    uPlantArray(3,(kk - BufferLength + 1):kk);    % FI-106 [sL/min]
                    ones(3,BufferLength); %dummy values --> in the actual rig, they will the pressure at injection point (PI105, PI106, PI107). Not used here
                    ones(3,BufferLength); %dummy values --> in the actual rig, they will the well temperature (TI101, TI102, TI103). Not used here
                    ones(3,BufferLength); %dummy values --> in the actual rig, they will be dP in a given pipe section (dPI101, dPI102, dPI103). Not used here
                    (measPlantArray(4,(kk - BufferLength + 2):kk + 1) - 1.01325)*10^3; % PI-101 [mbar g]
                    (measPlantArray(5,(kk - BufferLength + 2):kk + 1) - 1.01325)*10^3; % PI-102 [mbar g] 
                    (measPlantArray(6,(kk - BufferLength + 2):kk + 1) - 1.01325)*10^3; % PI-103 [mbar g]
                    uPlantArray(7,(kk - BufferLength + 1):kk) - 1.01325];     % PI-104 [bar g]
                
       
        % values of the input variables at the previus rig sampling time (dim = nu[7] X 1)
        P_vector = [uk(4); % CV101 opening [-]
                    uk(5); % CV102 opening [-]
                    uk(6); % CV103 opening [-]
                    1;%dummy values --> in the actual rig, they will the pump rotation. Not used here
                    uk(1);  % FI-104 [sL/min]
                    uk(2);  % FI-105 [sL/min]
                    uk(3)]; % FI-106 [sL/min]

        % values of the inputs (gas lift) of the last optimization run (dim = nQg[3] X 1)
        O_vector = uSPPlantArray(:,kk - nExec);
        
        % Run Labview/Matlab interface file
        LabViewMain

        flagArray = [flagArray, [SS;Estimation;Optimization]]; 
        ofArray = [ofArray, Result]; 
        thetaHatArray = [thetaHatArray, Parameter_Estimation']; 
        yEstArray = [yEstArray, State_Variables_Estimation']; 
        yOptArray = [yOptArray, State_Variables_Optimization']; 
        uOptArray = [uOptArray, Optimized_Air_Injection']; 

        uImpArray = [uImpArray, O_vector']; 
        
    else
         % update with dummy values
         flagArray = [flagArray, [0;0;0]];
         ofArray = [ofArray, 0];
         thetaHatArray = [thetaHatArray, zeros(1,6)'];
         yEstArray = [yEstArray, zeros(1,6)'];
         yOptArray = [yOptArray, zeros(1,6)'];
         uOptArray = [uOptArray, zeros(1,3)'];
         
         uImpArray = [uImpArray, zeros(1,3)'];
     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % updating input and parameter vectors %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1. Saving setpoints
    uSP = [O_vector(1);  %FIC-104 SP [sl/min]
           O_vector(2);  %FIC-105 SP [sl/min]
           O_vector(3)]; %FIC-106 SP [sl/min]
    uSPPlantArray = [uSPPlantArray, uSP];
   
    %2. Adding "controller" action
    % instead of coding a controller, we simply add a 5s delay + input noise
    if kk > 5
        uImple = uImple + (uSPPlantArray(:,kk - 5) - uImple)*exp(-parPlant.T/tau_C) + noise.input*randn(3,1);
    else
        uImple = uSPPlantArray(:,kk) + noise.input*randn(3,1);
    end
    
    %3. Saving actual implemented values
    uk = [uImple;
          disturbances.values(:,kk + 1)]; % for obtaining the full plant input vector, we need to add the values that correspond to the disturbances (feed)
     uPlantArray = [uPlantArray, uk];

     %4. Updating plant parameters according to pre-computed array
     %   the values are updated every 10s
     thetak = parProfile.thetaPlant(:,kk); 
     thetaPlantArray = [thetaPlantArray, thetak];
end

% save(name,'flagArray','ofArray','thetaHatArray','xEstArray','xOptArray','uOptArray','uImpArray'); 

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% plotting the data 
% checking sampling rate
markers = {'o','x','>'};
cc = {'b','k','g'};
leg = {'w_1','w_2','w_3'};
 
%% Inputs (Gas flowrate)
f = figure(1);
for well = 1:3
    subplot(3,1,well)
        plot(tgrid, uPlantArray(well,:),'b','Linewidth',1.5)
        hold on 
        plot(tgrid, uSPPlantArray(well,:),'k:','Linewidth',1.5)
        
        % chosen manually
        if well == 1
            ylim([2 3])
        elseif well == 2
            ylim([1 2.5])
        else
            ylim([2.5 4.5])
        end
        
        xticks(0:1:(1/60)*(nFinal - 1))
        xlim([0 (1/60)*(nFinal - 1)])

        xlabel('time [min]','FontSize',10)
        ylabel('Q_g [sL/min]','FontSize',10)
        
        name = ['Well ',num2str(well)];
        title(name,'FontSize',10)   
end
     
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!! If the parameters are set to zero, it means that the RTO   %
%    didn't run in this iteration                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for well = 1:3
    f = figure(2 + well);
    
    % Reservoir Valve Parameters
    subplot(4,1,1)
        hold on 
        plot(tgrid(2:end), yEstArray(well,:),'bo','linestyle','none','markersize',5)
        plot(tgrid, measPlantArray(well,:),'k:','Linewidth',1.5)
       
        xticks(0:5:(1/60)*(nFinal - 1))
        xlim([0 (1/60)*(nFinal - 1)])
        % chosen manually
        if well == 1
            ylim([5 10])
        elseif well == 2
            ylim([3 8])
        else
            ylim([5 10])
        end
    
        xlabel('time [min]','FontSize',10)
        ylabel('[L/min]','FontSize',10)
        title('Liquid Flowrate','FontSize',10)
        
    subplot(4,1,2)
        hold on
        plot(tgrid(2:end), thetaHatArray(well,:),'bo','linestyle','none','markersize',5)
        plot(tgrid, thetaPlantArray(well,:),'k:','Linewidth',1.5)
        
        % chosen manually
        if well == 1
            ylim([0.01 0.5])
        elseif well == 2
            ylim([0.01 0.2])
        else
            ylim([0.01 0.5])
        end
        
        xticks(0:5:(1/60)*(nFinal - 1))
        xlim([0 (1/60)*(nFinal - 1)])

        xlabel('time [min]','FontSize',10)
        title('Reservoir Parameters','FontSize',10)

    subplot(4,1,3)
        hold on
        plot(tgrid(2:end), yEstArray(3 + well,:),'bo','linestyle','none','markersize',5)
        plot(tgrid, measPlantArray(3 + well,:),'k:','Linewidth',1.5)
        
        xticks(0:5:(1/60)*(nFinal - 1))
        xlim([0 (1/60)*(nFinal - 1)])
        ylim([0.9 1.1])

        xlabel('time [min]','FontSize',10)
        ylabel('[mbar G]','FontSize',10)
        title('Top pressure','FontSize',10)
        
   subplot(4,1,4)
        hold on
        plot(tgrid(2:end), thetaHatArray(3 + well,:),'bo','linestyle','none','markersize',5)
        plot(tgrid, thetaPlantArray(3 + well,:),'k:','Linewidth',1.5)  
        
        legend({'Estimated','True'},'Position',[0.727767857142857 0.214091213084366 0.167261904761905 0.0710317460317459]);
        
        % chosen manually
        if well == 1
            ylim([0.4 1])
        elseif well == 2
            ylim([0.4 1])
        else
            ylim([0.5 1.1])
        end
        
        xticks(0:5:(1/60)*(nFinal - 1))
        xlim([0 (1/60)*(nFinal - 1)])

        xlabel('time [min]','FontSize',10)
        title('Valve Parameters','FontSize',10)

end

%% SS Detection
f = figure(6);

subplot(3,1,1)
    hold on
    plot(tgrid,uPlantArray(4,:),'k:','Linewidth',1.5)
    plot(tgrid,uPlantArray(5,:),'k:','Linewidth',1.5)
    plot(tgrid,uPlantArray(6,:),'k:','Linewidth',1.5)
    
    hold off

    xticks(0:1:(1/60)*(nFinal - 1))
    xlim([0 (1/60)*(nFinal - 1)])
    
    legend({'CV101','CV102','CV103'},'Position',[0.77 0.82 0.135 0.10]);

    title('Disturbances')
    ylabel('valve open. [%]')
    xlabel('time [min]')

        
% RTO executions    
subplot(3,1,3)   
    plot(tgrid(2:end),flagArray(1,:),'x','MarkerSize',5) 
        
    xticks(0:1:(1/60)*(nFinal - 1))
    xlim([0 (1/60)*(nFinal - 1)])
    ylim([-0.5 1.5])
    yticks([0 1])
    yticklabels({'No','Yes'})
    
    title('RTO executions')
      
    
% Liquid measurements    
subplot(3,1,2)
    hold on
    for ii = 1:3
       plot(tgrid, measPlantArray(ii,:),'Color',cc{ii},'Linewidth',1.5)    
    end
    
    hold off
    legend({'Well 1','Well 2','Well 3'},'Location','northwest')
    
    xlim([0 (1/60)*(nFinal - 1)])
    xticks(0:1:(1/60)*(nFinal - 1))

    ylabel('Q_l [L/min]')
    xlabel('time [min]')
