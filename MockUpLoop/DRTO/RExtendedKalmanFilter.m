function [dxHatkk,zHatkk,pHatkk,Pkk,qThres] = RExtendedKalmanFilter(yk,uk_1,dxk_1k_1,zk_1k_1,p_k_1k_1,PNk_1k_1,F,S_xx,S_zz,S_xz,S_xp,S_zp,ekf,par)
%    Uses an reduced extended Kalman filter for estimating the model states and parameters
%    We set the parameters qThres and ekf.lambdaThres such that the rEKF
%    calculations become equivalent to the EKF

% Inputs:
%    yk, uk_1 = measurements + system inputs
%    dxk_1k_1,zk_1k_1,p_k_1k_1 = previously estimated states and parameters
%    PNk_1k_1 = previously computed estimate covariance matrix
%    F = dynamic model integrator
%    S_xx,S_zz,S_xz,S_xp,S_zp = sensitivty equations
%    ekf = ekf configuration parameters
%    par = system parameters

% Outputs:
%   pHatkk = estimated parameters
%   dxHatkk, zHatkk = estimated states
%   Pkk = updated value of the estimate covariance
%   qThres  = number of parameters + states estimated --> not used here

% Other m-files required: none
% MAT-files required: none

% addpath ('C:\Users\lab\Documents\casadi-windows-matlabR2016a-v3.4.5')        
import casadi.*

% ===================================
%     Integrating results
% ===================================
    % we use the nonlinear model to evolve the states
    Fend = F('x0',dxk_1k_1,'z0',zk_1k_1,'p',[uk_1;p_k_1k_1;par.T]);
    %extracting solution
    dxkk_1 = full(Fend.xf);
    zkk_1 = full(Fend.zf);
    
% ===================================
%     Sensitivitits
% ===================================
    %states
    F_xx = full(S_xx(dxk_1k_1,zk_1k_1,[uk_1;p_k_1k_1;par.T]));
    F_zz = full(S_zz(dxk_1k_1,zk_1k_1,[uk_1;p_k_1k_1;par.T]));
    F_xz = full(S_xz(dxk_1k_1,zk_1k_1,[uk_1;p_k_1k_1;par.T]));
    F_zx = F_xz';    
    
    FSk = [F_xx, F_xz; F_zx, F_zz];
    HSk = [zeros(length(yk),length(dxk_1k_1)),par.H];

    %inputs
    F_xp = full(S_xp(dxk_1k_1,zk_1k_1,[uk_1;p_k_1k_1;par.T]));
    F_zp = full(S_zp(dxk_1k_1,zk_1k_1,[uk_1;p_k_1k_1;par.T]));
    FPk = [F_xp;F_zp];
    FPk(:,1:length(uk_1)) = [];%excluding derivative in relation to feed inputs
    FPk(:,end) = [];%excluding derivative in relation to par.T
    
% ===================================
%     Filter
% ===================================
    %Arranging filter matrices - extended state vector
    Fk_1 = [FSk, FPk; zeros(length(p_k_1k_1),(length(dxk_1k_1) + length(zkk_1))), eye(length(p_k_1k_1))]; 
    
    %output transition
    Hk = [HSk, zeros(length(yk),length(p_k_1k_1))]; 
 
    %preparing extended vector
    xkk_1 = [dxkk_1; zkk_1; p_k_1k_1];

    %%%%%%%%%%%%%%%%%%%%
    % Filter Equations %
    %%%%%%%%%%%%%%%%%%%%
    %predicted covariance estimate
    Pkk_1 = Fk_1*PNk_1k_1*Fk_1'+ ekf.Qe;

    %%
    % N.B.: if EKF is used instead of rEKF, this section is not relevant
    %Reduced Filtering - spectral deconposition
    [V,D] = eig(Pkk_1);

    %main diagonal - abs value of eigenvalues
    d = diag(D);

    %sorting eigenvalues - descend
    %[e,i] = sort(d,'descend');
    [e,i] = sort(d);

    %initializing the threshold. The initial value is equal to the total number
    %of states + parameters
    qThres = size(ekf.Qe,2);
    
    %initializing counter
    l = 1;
    while l <= length(e)
        if e(l) > ekf.lambdaThres
            %checking if the eigenvalue is greater than the threshold. If it
            %is, the number of directions used changes
            qThres = l;
            l = length(e);%exiting the loop
        end
        l = l + 1;
    end

    %obtaining reduced covariance matrix
    temp = zeros(size(Pkk_1,1));
    for j = 1:qThres
        temp = temp + e(j)*V(:,i(j))*V(:,i(j))';
    end
    Pkk_1r = temp;

    temp = [];

    %%
    %Updating
    %Innovation covariance
    Sk = Hk*Pkk_1r*Hk' + ekf.R;

    %Kalman Gain
    Kk = Pkk_1r*Hk'*pinv(Sk);

    %Updating covariance estimate
    %Galo's upload system! Reduced Matrix is used only in the gain calculations 
    Pkk = (eye(length(dxk_1k_1) + length(zk_1k_1) + length(p_k_1k_1)) - Kk*Hk)*Pkk_1*((eye(length(dxk_1k_1) + length(zk_1k_1) + length(p_k_1k_1)) - Kk*Hk))' + Kk*ekf.R*Kk';

    %guaranteeing symmetric matrix
    temp = Pkk;

    for i = 1:size(temp,1)
        for j = (i + 1):size(temp,2)
            temp(j,i) = Pkk(i,j);
        end
    end
    Pkk = temp;

    %Measurement residual
    yHat = Hk*xkk_1;

    %Update state estimate
    xHatkk = xkk_1 + Kk*(yk - yHat); %state variables are normalized
    
    % ===================================
    %     Results
    % ===================================
    %dividing extended vector
    dxHatkk = xHatkk(1:length(dxk_1k_1));%differential
    zHatkk = xHatkk(length(dxk_1k_1) + 1:length(dxk_1k_1) + length(zk_1k_1));%algebraic
    pHatkk = xHatkk(length(dxk_1k_1) + length(zk_1k_1) + 1:end);%parameters

end

