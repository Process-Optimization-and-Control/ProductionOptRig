  function [flagSystem,pArray] = SSDetection(yMed)
  %    Checks data window xMed and detects steady-state using linear method
  
  % Inputs:
  %    yMed = representative measurements for SSD
  %
  % Outputs:
  %   flagSystem: flac containing the SS flag (== 0 not at SS, == 1 at SS)
  %   pArray: pValue of the linear coefficients of the linear model
  
  % Other m-files required: none
  % MAT-files required: none
  
    % normalize measurements
    yMed_norm = normalize(yMed);
    
    % create time array for linear model:
    % y = A*t + b --> x_norm = A*t
    time = 1:size(yMed,2);
  
    % arrays for the case where there are more than 1 measurements in yMed
    flagArray = [];
    pArray = [];
    
    for ii = 1:size(yMed_norm,1) %checking all measurements
                                 % all of them need to be at SS!
        %creating a table 
        y_norm = yMed_norm(ii,:);
        tbl = table(y_norm',time');
        
        %estimating linear model
        mdl = fitlm(tbl,'Var1 ~ Var2');
        
        %The model display includes the p-value for the t-statistic for each coefficient to test the null hypothesis that the corresponding coefficient is zero.
        [p_time,~,~] = coefTest(mdl,[0 1]);
        pArray = [pArray; p_time];
        
        if p_time < 0.1 % significancy of the test
            %Failed! Most likely not at steady-state
            flagArray = [flagArray; 0];
        else
            %Passed
            flagArray = [flagArray; 1];
        end
    end
    
    %product of array entries - complete system
    flagSystem = prod(flagArray);
    
    
  end
