clc
clear 
close all

% Number of runs
nR = 2;

%load files
rd{1,1} = load('SSRTO_results_trial_1');
rd{1,2} = load('SSRTO_results_trial_2');
rd{1,3} = load('SSRTO_results_trial_3');

rd{2,1} = load('HRTO_results_trial_1');
rd{2,2} = load('HRTO_results_trial_2');

rd{3,1} = load('DRTO_results_trial_2');
rd{3,2} = load('DRTO_results_trial_3');

% Optimization sampling time
sT = 10;

% to control figure printing
% false - do not print
% true - print
pPDF = true;
pTIFF = false;
pr = [pPDF, pTIFF];

mark = {'-','-','-'}; % well 1 | 2 |3
mark2 = {'o','none','none'}; % well 1 | 2 |3
line2 = {'none','-','-'}; % well 1 | 2 |3

%\definecolor{HRTO}{RGB}{14,76,146} % YALE #0E4C92
%\definecolor{DRTO}{RGB}{255,36,0} % SCARLET #FF2400 
%\definecolor{SSRTO}{RGB}{207,152,18} % TANGERINE #CF9812A
colorMethod = 1/255*[207,152,18;
                     14,76,146;
                     255,36,0];
                 
colorMap = [0.90,0.90,0.90; %light gray
            0.95,0.95,0.95]; % gelo
cmap = repmat(colorMap,5,1);

colorWell = {'k','b','r'}; % well 1 | 2 |3
inputTitle = {'Well 1','Well 2','Well 3'};
methodTitle = {'SSRTO','ROPA','DRTO'};
parameterTitle_1 = {'\theta_{res,1}','\theta_{res,2}','\theta_{res,3}'};
parameterTitle_2 = {'\theta_{val,1}','\theta_{val,2}','\theta_{val,3}'};

%%
%%%%%%%% 
% DATA %
%%%%%%%%
% final order of the data is 
% DATA{1} (every second)
 % 1: Relative time [s]	
 % 2: FIC-104 setpt [sl/min]	
 % 3: FIC-104 [sl/min]	
 % 4: FIC-105 setpt [sl/min]	
 % 5: FIC-105 [sl/min]	
 % 6: FIC-106 setpt [sl/min]	
 % 7: FIC-106 [sl/min]	
 % 8: CV-101 setpt [l/min]	
 % 9: FI-101 [l/min]	
 % 10: CV-102 setpt [l/min]	
 % 11: FI-102 [l/min]
 % 12: CV-103 setpt [l/min]	
 % 13: FI-103 [l/min]	
 % 14: dP-101 [mbar D]	
 % 15: dP-102 [mbar D]	
 % 16: dP-103 [mbar D]	
 % 17: CV-107 setpt [mbar G]	
 % 18: PI-101 [mbar G]	
 % 19: CV-108 setpt [mbar G]	
 % 20: PI-102 [mbar G]	
 % 21: CV-109 setpt [mbar G]	
 % 22: PI-103 [mbar G]	
 % 23: Pump output pressure setpt [bar g]	
 % 24: PI-104 [bar G]
 % 25: TI-101 [oC]
 % 26: TI-102 [oC]
 % 27: TI-103 [oC]
 % 28: CV-101 current [A]	
 % 29: CV-102 current [A]	
 % 30: CV-103 current [A]	
 % 31: CV-107 current [A]	
 % 32: CV-108 current [A]	
 % 33: CV-109 current [A]	
 % 34: Pump ctrl current [A]
 
% DATA{2} (every minute)
 % 1: SS 
 % 2: EstimationError 
 % 3: OptmizationError 
 % 4: OF 
 % 5: Theta-w1 
 % 6: Theta-w2 
 % 7: Theta-w3 
 % 8: Valve-w1 
 % 9: Valve-w2 
 %10: Valve-w3 
 %11: WroEstimated-w1 
 %12: WroEstimated-w2 
 %13: WroEstimated-w3 
 %14: PrhEstimated-w1 
 %15: PrhEstimated-w2 
 %16: PrhEstimated-w3 
 %17: WroOptimzed-w1 
 %18: WroOptimzed-w2 
 %19: WroOptimzed-w3 
 %20: PrhOptimzed-w1 
 %21: PrhOptimzed-w2 
 %22: PrhOptimzed-w3 
 %23: u0-w1 
 %24: u0-w2 
 %25: u0-w3
 
 %% mean and std calculations 
 % three independent runs with each method. Here, we take the mean and
 % average - the times are not exactly the same, but we approximate it as
 % being the same 

 %dummy
 minRun = zeros(1,2);
 minLength = 1e5; 

 %checking the minimum time length of all the runs
 for ii = 1:3 % methods
     for jj = 1:nR % runs
         tempLength = length(rd{ii,jj}.DATA{1}.B(1,:));
        if  tempLength < minLength
            minLength = tempLength;
            minRun(1) = ii;
            minRun(2) = jj;
        end
         
     end
 end
 
% computing the average and std of all the measurements in the three runs

% generating empty variable for saving the data
for ii = 1:3
    mLV{ii}.mean = zeros(34,minLength);
    mLV{ii}.std = zeros(34,minLength);
    mMA{ii}.mean = zeros(25,floor(minLength/sT));
    mMA{ii}.std = zeros(25,floor(minLength/sT));
    devMeanProfitLV{ii} = [];
end

% for LABVIEW DATA
time = 1:minLength;



for kk = 1:length(time)
   for ii = 1:3 %methods
      temp = [];
      for jj = 1:nR         
          temp = [temp, rd{ii,jj}.DATA{1}.B(:,time(kk))];
      end
                 
      mLV{ii}.mean(:,kk) = mean(temp,2);
      mLV{ii}.std(:,kk) = std(temp,[],2);
      
      temp2 = [];
      for jj = 1:nR         
          temp2 = [temp2, 20*rd{ii,jj}.DATA{1}.B(9,time(kk)) + 10*rd{ii,jj}.DATA{1}.B(11,time(kk)) + 30*rd{ii,jj}.DATA{1}.B(13,time(kk))];
      end
      
      mLV{ii}.profit(:,kk) = mean(temp2,2);
      mLV{ii}.profit_std(:,kk) = std(temp2,[],2);
      
      devMeanProfitLV{ii} = [devMeanProfitLV{ii}, (temp2 - mLV{ii}.profit(:,kk)).^2];
   end
   
   
end

% for MATLAB DATA
timeRTO = 1:sT:minLength;
for kk = 1:length(timeRTO)
   for ii = 1:3
      temp = [];
      for jj = 1:nR         
          temp = [temp, rd{ii,jj}.DATA{2}.B(:,kk)];
      end
                
       mMA{ii}.mean(:,kk) = mean(temp,2);
       mMA{ii}.std(:,kk) = std(temp,[],2);
       
     
      temp2 = [];
      for jj = 1:nR         
          temp2 = [temp2, 20*rd{ii,jj}.DATA{2}.B(17,kk) + 10*rd{ii,jj}.DATA{2}.B(18,kk) + 30*rd{ii,jj}.DATA{2}.B(19,kk)];
      end
      
      mMA{ii}.profit(:,kk) = mean(temp2,2);
      mMA{ii}.profit_std(:,kk) = std(temp2,[],2);
       
   end
end

% No RTO
noRTO = load('no_RTO');
noRTO_profit = 20*noRTO.DATA{1}.B(9,:) + 10*noRTO.DATA{1}.B(11,:) + 30*noRTO.DATA{1}.B(13,:);
ntNoRTO = length(noRTO.DATA{1}.time);

proDifPerArray = [];
meanProfDif= [];

for ii = 1:3
    temp = 100*(mLV{ii}.profit - noRTO_profit(1:minLength))./abs(noRTO_profit(1:minLength));
    
    % applying moving average
    proDifPerArray = [proDifPerArray; medfilt1(temp,60)];
    
    meanProfDif= [meanProfDif, mean(temp)];
end

% % meanProfDif
% temp1 = cumsum(proDifPerArray(1,:)/60);
% temp2 = cumsum(proDifPerArray(2,:)/60);
% temp3 = cumsum(proDifPerArray(3,:)/60);
% 
% temp1(end)
% temp2(end)
% temp3(end)
% 
% temp4 = mean(proDifPerArray(1,:))
% temp5 = mean(proDifPerArray(2,:))
% temp6 = mean(proDifPerArray(3,:))

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
%% Disturbances vs. Inputs
f = figure(1);

for well = 1:3
    
subplot(3,1,well)    
    hold on
    
    %plotting disturbance regions
    for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
        X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
        Y = ones(1,length(X))*120;
        h = area(X,Y,'LineStyle','none','HandleVisibility','off');
        h(1).FaceColor = cmap(kk,:);
    end

        for ii = 1:3 %methods mean
            plot(time/60,mLV{ii}.mean(3 + 2*(well - 1),:),'Color',colorMethod(ii,:),'LineStyle',mark{ii},'Linewidth',1.5)   
        end
        
        yline(1,':','LineWidth',1,'HandleVisibility','off');
        yline(5.5,':','LineWidth',1,'HandleVisibility','off');
        hold off

        % limits
        xlim([0 time(end)/60])
        
        % chosen manually
        if well == 1
            ylim([2 3])
        elseif well == 2
            ylim([1 2.5])
        else
            ylim([2.5 4.5])
        end
        
        if well == 2
            legend(methodTitle,'Position',[0.75 0.50 0.16 0.12]);
        end
        xlabel('time [min]')
        ylabel('Q_{g} [sL/min]')
        
        tit = ['Well ',num2str(well)];
        title(tit)
             
end

if pr(1)
    save_name = 'Results_Inputs.pdf';
    print(f,save_name,'-dpdf')
end
if pr(2)
    save_name = 'Results_Inputs_.tif';
    print(f,'-r1200','-dtiff',save_name);
end

%% Inputs variability
for well = 1:3
    f = figure(well + 3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Different implementations 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:3 %methods
        subplot(3,1,ii)
        hold on
        
        %plotting disturbance regions
        for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
            X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
            Y = ones(1,length(X))*120;
            h = area(X,Y,'LineStyle','none','HandleVisibility','off');
            h(1).FaceColor = cmap(kk,:);
        end
        
        % mean
        plot(time/60,mLV{ii}.mean(3 + 2*(well - 1),:),'Color',colorMethod(ii,:),'LineStyle',mark{ii},'Linewidth',1.5)
        % plus/minus std
        plot(time/60,mLV{ii}.mean(3 + 2*(well - 1),:) + 2*mLV{ii}.std(3 + 2*(well - 1),:),'Color',colorMethod(ii,:),'LineStyle',':','Linewidth',1)
        plot(time/60,mLV{ii}.mean(3 + 2*(well - 1),:) - 2*mLV{ii}.std(3 + 2*(well - 1),:),'Color',colorMethod(ii,:),'LineStyle',':','Linewidth',1)

        hold off
        
        % limits
        xlim([0 time(end)/60])
        
        % chosen manually
        if well == 1
            ylim([2 3.5])
        elseif well == 2
            ylim([1 2.5])
        else
            ylim([2.5 4.5])
        end
        
        legend({'\mu','\mu + 2\sigma','\mu - 2\sigma'},'Location','best');

        tit = ['Well ',num2str(well),': ',methodTitle{ii}];
        title(tit);
        xlabel('time [min]')
        ylabel('Q_{g} [sL/min]')
    end
      
    if pr(1)
        save_name = ['Results_InputsVar_Well_',num2str(well),'.pdf'];
        print(f,save_name,'-dpdf')
    end
    if pr(2)
        save_name = ['Results_InputsVar_Well_',num2str(well),'.tif'];
        print(f,'-r1200','-dtiff',save_name);
    end
           
end

%% Profit variability
f = figure(7);

for ii = 1:3
    
    subplot(3,1,ii)
    hold on
    
    %plotting disturbance regions
    for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
        X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
        Y = ones(1,length(X))*600;
        h = area(X,Y,'LineStyle','none','HandleVisibility','off');
        h(1).FaceColor = cmap(kk,:);
    end
    
    % mean
    plot(time/60,mLV{ii}.profit,'Color',colorMethod(ii,:),'LineStyle',mark{ii},'Linewidth',1.5)
    % plus/minus std
    plot(time/60,mLV{ii}.profit + 2*sqrt(mean(devMeanProfitLV{ii})),'Color',colorMethod(ii,:),'LineStyle',':','Linewidth',0.75)
    plot(time/60,mLV{ii}.profit - 2*sqrt(mean(devMeanProfitLV{ii})),'Color',colorMethod(ii,:),'LineStyle',':','Linewidth',0.75)
    
    hold off
    
    % limits
    xlim([0 time(end)/60])
    ylim([280 520])
    
    legend({'\mu','\mu + 2\sigma','\mu - 2\sigma'},'Location','best');
    
    title(methodTitle{ii});
    xlabel('time [min]')
    ylabel('J_{inst.} [$/min]')
end
      
    if pr(1)
        save_name = 'Results_ProfitVar.pdf';
        print(f,save_name,'-dpdf')
    end
    if pr(2)
        save_name = 'Results_ProfitVar.tif';
        print(f,'-r1200','-dtiff',save_name);
    end


%% SS Detection
f = figure(8);

subplot(3,1,1)
    hold on
    
    %plotting disturbance regions
    for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
        X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
        Y = ones(1,length(X))*120;
        h = area(X,Y,'LineStyle','none','HandleVisibility','off');
        h(1).FaceColor = cmap(kk,:);
    end

    plot(rd{minRun(1),minRun(2)}.DATA{1}.time/60,100*(rd{minRun(1),minRun(2)}.DATA{1}.B(28,:) - 0.004)./(0.02 - 0.004),'k:','Linewidth',1.5)
    plot(rd{minRun(1),minRun(2)}.DATA{1}.time/60,100*(rd{minRun(1),minRun(2)}.DATA{1}.B(29,:) - 0.004)./(0.02 - 0.004),'k--','Linewidth',1.5)
    plot(rd{minRun(1),minRun(2)}.DATA{1}.time/60,100*(rd{minRun(1),minRun(2)}.DATA{1}.B(30,:) - 0.004)./(0.02 - 0.004),'k-.','Linewidth',1.5)
    
    hold off

    xlim([0 rd{minRun(1),minRun(2)}.DATA{1}.time(end)/60])
    ylim([10 90])
    legend({'CV101','CV102','CV103'},'Position',[0.77 0.82 0.135 0.10]);

    title('Disturbances')
    ylabel('valve open. [%]')
    xlabel('time [min]')

        
% subplot(3,1,2)
%     hold on 
%     
%     %plotting disturbance regions
%     for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
%         X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
%         Y = ones(1,length(X))*120;
%         h = area(X,Y,'LineStyle','none','HandleVisibility','off');
%         h(1).FaceColor = cmap(kk,:);
%     end
%     
%     plot(rd{1,1}.DATA{1}.time/60,3*rd{1,1}.DATA{1}.B(9,:),'r:','Linewidth',1.5)
%     plot(rd{1,1}.DATA{1}.time/60,2*rd{1,1}.DATA{1}.B(11,:),'k--','Linewidth',1.5)
%     plot(rd{1,1}.DATA{1}.time/60,1*rd{1,1}.DATA{1}.B(13,:),'b-.','Linewidth',1.5)
% 
%     
%     hold off
%     
%     xlim([0 time(end)/60])
%     ylim([10 14])
%     xlabel('time [min]')
%     ylabel('Q_{l} [L/min]')
% 
%     title('Representative Measurements (liquid flowrates)')
%     
subplot(3,1,2)
    hold on 
    plot(timeRTO/60,1*rd{1,1}.DATA{2}.B(1,:),'x','MarkerSize',5)
    plot(timeRTO/60,2*rd{1,2}.DATA{2}.B(1,:),'x','MarkerSize',5)
    plot(timeRTO/60,3*rd{1,3}.DATA{2}.B(1,:),'x','MarkerSize',5)
    grid on 
    box on
    
    xlim([0 timeRTO(end)/60])
    ylim([0.5 3.5])
    yticks(1:3)
    yticklabels({'Run 1','Run 2','Run 3'})
    
    title('Steady-state detection flag')
    xlabel('time [min]')

    if pr(1)
        save_name = 'Results_SSDet.pdf';
        print(f,save_name,'-dpdf')
    end
    if pr(2)
        save_name = 'Results_SSDet.tif';
        print(f,'-r1200','-dtiff',save_name);
    end

%% Disturbances vs. Parameter estimates (Reservoir Valves)
 f = figure(9);

for well = 1:3
subplot(3,1,well)    
    hold on
    
    %plotting disturbance regions
    for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
        X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
        Y = ones(1,length(X))*120;
        h = area(X,Y,'LineStyle','none','HandleVisibility','off');
        h(1).FaceColor = cmap(kk,:);
    end

        for ii = 1:3 %methods mean
            isNZ=(~rd{ii,1}.DATA{2}.B((4 + well),:)==0);           % addressing logical array of nonzero elements
            plot(timeRTO(isNZ)/60,rd{ii,1}.DATA{2}.B((4 + well),isNZ),...
                'Color',colorMethod(ii,:),...
                'Linewidth',1.5,...
                'Marker',mark2{ii},...
                'LineStyle',line2{ii})
            
            
%             isNZ=(~mMA{ii}.mean((4 + well),:)==0);           % addressing logical array of nonzero elements
%             plot(timeRTO(isNZ)/60,mMA{ii}.mean((4 + well),isNZ),...
%                 'Color',colorMethod(ii,:),...
%                 'Linewidth',1.5,...
%                 'Marker',mark2{ii},...
%                 'LineStyle',line2{ii})

        end
        
        %yline(1,':','LineWidth',1,'HandleVisibility','off');
        %yline(5.5,':','LineWidth',1,'HandleVisibility','off');
        hold off

        % limits
        xlim([0 time(end)/60])
        
        % chosen manually
        if well == 1
            ylim([0 0.5])
        elseif well == 2
            ylim([0 0.2])
        else
            ylim([0 0.5])
        end
        
        if well == 2
            legend(methodTitle,'Position',[0.75 0.82 0.16 0.12]);
        end
        
        xlabel('time [min]')
        ylabel('\theta_{res}')
        
        tit = ['Well ',num2str(well)];
        title(tit)
              
end

if pr(1)
    save_name = 'Results_ParReservoir.pdf';
    print(f,save_name,'-dpdf')
end
if pr(2)
    save_name = 'Results_ParReservoir_.tif';
    print(f,'-r1200','-dtiff',save_name);
end

%% Disturbances vs. Parameter estimates (Top Valves)
f = figure(12);

for well = 1:3
subplot(3,1,well)   
    hold on
    
    %plotting disturbance regions
    for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
        X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
        Y = ones(1,length(X))*120;
        h = area(X,Y,'LineStyle','none','HandleVisibility','off');
        h(1).FaceColor = cmap(kk,:);
    end

        for ii = 1:3 %methods mean
            isNZ=(~rd{ii,2}.DATA{2}.B((7 + well),:)==0);           % addressing logical array of nonzero elements
            plot(timeRTO(isNZ)/60,rd{ii,2}.DATA{2}.B((7 + well),isNZ),...
                    'Color',colorMethod(ii,:),...
                    'Linewidth',1.5,...
                    'Marker',mark2{ii},...
                    'LineStyle',line2{ii})                           

                
%             isNZ=(~mMA{ii}.mean((4 + well),:)==0);           % addressing logical array of nonzero elements
%             plot(timeRTO(isNZ)/60,mMA{ii}.mean((4 + well),isNZ),...
%                     'Color',colorMethod(ii,:),...
%                     'Linewidth',1.5,...
%                     'Marker',mark2{ii},...
%                     'LineStyle',line2{ii}) 
        end
        
        %yline(1,':','LineWidth',1,'HandleVisibility','off');
        %yline(5.5,':','LineWidth',1,'HandleVisibility','off');
        hold off

        % limits
        xlim([0 time(end)/60])
        
        % chosen manually
        if well == 1
            ylim([0.4 1])
        elseif well == 2
            ylim([0.4 1])
        else
            ylim([0.5 1.1])
        end
        
%         if well == 2
%             legend(methodTitle,'Position',[0.75 0.50 0.16 0.12]);
%         end
        
        xlabel('time [min]')
        ylabel('\theta_{top}')
        
        tit = ['Well ',num2str(well)];
        title(tit)
                
end

if pr(1)
    save_name = 'Results_ParTop.pdf';
    print(f,save_name,'-dpdf')
end
if pr(2)
    save_name = 'Results_ParTop.tif';
    print(f,'-r1200','-dtiff',save_name);
end

%% Disturbances vs. Delta Input
f = figure(15);
% subplot(4,1,1)
%     hold on
%     
%     %plotting disturbance regions
%     for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
%         X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
%         Y = ones(1,length(X))*120;
%         h = area(X,Y,'LineStyle','none','HandleVisibility','off');
%         h(1).FaceColor = cmap(kk,:);
%     end
% 
%     plot(rd{minRun(1),minRun(2)}.DATA{1}.time/60,100*(rd{minRun(1),minRun(2)}.DATA{1}.B(28,:) - 0.004)./(0.02 - 0.004),'k:','Linewidth',1.5)
%     plot(rd{minRun(1),minRun(2)}.DATA{1}.time/60,100*(rd{minRun(1),minRun(2)}.DATA{1}.B(29,:) - 0.004)./(0.02 - 0.004),'k--','Linewidth',1.5)
%     plot(rd{minRun(1),minRun(2)}.DATA{1}.time/60,100*(rd{minRun(1),minRun(2)}.DATA{1}.B(30,:) - 0.004)./(0.02 - 0.004),'k-.','Linewidth',1.5)
%     
%     hold off
% 
%     xlim([0 rd{minRun(1),minRun(2)}.DATA{1}.time(end)/60])
%     ylim([10 90])
%     legend({'well 1','well 2','well 3'},'Position',[0.792261904761905 0.840429079345783 0.135119047619048 0.102380952380952]);
% 
%     title('Disturbances')
%     ylabel('valve open. [%]')
% 
% for well = 1:3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % Different wells 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(4,1,1 + well)    
%     hold on
%     
%     %plotting disturbance regions
%     for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
%         X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
%         Y = ones(1,length(X))*120;
%         h = area(X,Y,'LineStyle','none','HandleVisibility','off');
%         h(1).FaceColor = cmap(kk,:);
%     end
% 
%         for ii = 1:3 %methods mean
%             
%             nData = size(mLV{ii}.mean,2);
%             deltaU = mLV{ii}.mean(2 + 2*(well - 1),2:nData) - mLV{ii}.mean(2 + 2*(well - 1),1:nData-1);
%             timeDu = 2:nData;
%             
%             stairs(timeDu/60,deltaU,...
%                     'Color',colorMethod(ii,:),...
%                     'Linewidth',1.5)                           
%         end
%         
%         %yline(1,':','LineWidth',1,'HandleVisibility','off');
%         %yline(5.5,':','LineWidth',1,'HandleVisibility','off');
%         hold off
% 
%         % limits
%         xlim([0 timeDu(end)/60])
%         
%         % chosen manually
%         ylim([0 0.2])
% 
%         %legend(methodTitle,'Location','best');
%         xlabel('time [min]')
%         ylabel('\Delta u [sL/min]')
%         
%         tit = ['Well ',num2str(well)];
%         title(tit)
%                   
% end

for well = 1:3

    for ii = 1:3
        nData = size(mLV{ii}.mean,2);
        deltaU = mLV{ii}.mean(2 + 2*(well - 1),2:nData) - mLV{ii}.mean(2 + 2*(well - 1),1:nData-1);
        
        subplot(1,3,well), 
            hold on
            plot(ii*ones(1,length(deltaU)),deltaU,'Color',colorMethod(ii,:),'marker','+','linestyle','none','MarkerSize',5)
        
    end
    box on 
    grid on
    
    ylim([-0.7 0.7])
    xlim([0.5 3.5])
    xticks(1:3)
    xticklabels(methodTitle)
    xtickangle(45)
    
    ylabel('\Delta Q_{gl} [sL/min]')
    title(inputTitle{well})

end
hold off

if pr(1)
    save_name = 'Results_DeltaU.pdf';
    print(f,save_name,'-dpdf')
end
if pr(2)
    save_name = 'Results_DeltaU.tif';
    print(f,'-r1200','-dtiff',save_name);
end

%% Disturbances vs. profit difference (naive)
f = figure(16);

for ii = 1:3
    subplot(4,1,ii)
    hold on
    %plotting disturbance regions
    for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
        X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
        Y = ones(1,length(X))*10;
        h = area(X,Y,'LineStyle','none','HandleVisibility','off');
        h(1).FaceColor = cmap(kk,:);
        h.BaseValue = - 5;
    end
    
    plot(time/60,proDifPerArray(ii,:),'Color',colorMethod(ii,:),'LineWidth',1.5)
    
    yline(0,'k:','LineWidth',1.0);
    %     yline(5,'k:','LineWidth',1.0,'HandleVisibility','off');
    %     yline(-5,'k:','LineWidth',1.0,'HandleVisibility','off');
    
    grid on
    
    %legend({'100(\psi - \psi_{naive})/\psi_{naive}','Filtered value','Reference'},'Location','northwest')
    title(methodTitle(ii))
    xlim([0 time(end)/60])
    ylim([-5 5])
    xlabel('time [min]')
    ylabel('J_{dif} [% /min]')

    
end

subplot(4,1,4)
hold on
%plotting disturbance regions
for kk = 1:length(rd{minRun(1),minRun(2)}.distArray) - 1
    X = [rd{minRun(1),minRun(2)}.distArray(kk)/60, rd{minRun(1),minRun(2)}.distArray(kk + 1)/60];
    Y = ones(1,length(X))*3000;
    h = area(X,Y,'LineStyle','none','HandleVisibility','off');
    h(1).FaceColor = cmap(kk,:);
end

for ii = 1:3
    %divide by 60 to adjust time scale
    plot(time/60,cumsum(proDifPerArray(ii,:)/60),'Color',colorMethod(ii,:),'LineWidth',1.5)
    
%     yline(0,'k:','LineWidth',1.0);
%     %     yline(5,'k:','LineWidth',1.0,'HandleVisibility','off');
%     %     yline(-5,'k:','LineWidth',1.0,'HandleVisibility','off');
%         
    %legend({'100(\psi - \psi_{naive})/\psi_{naive}','Filtered value','Reference'},'Location','northwest')
    title('Cummulative value')
    xlim([0 time(end)/60])
    ylim([0 40])
    yticks(0:10:40)
    xlabel('time [min]')
    ylabel('\Sigma J_{dif} [%]')

    
end

hold off

if pr(1)
    print(f,'Results_Profit_Dif','-dpdf')
end
if pr(2)
    print(f,'-r1200','-dtiff','Results_Profit_Dif.tif');
end


