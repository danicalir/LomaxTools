function [Ase,Bse,xrse,A_se,B_se,xr_se] = lomaxbootstrap(A_opt,B_opt,theta,CensorPoint,Precision,varargin)

% Calculates bootstrapped standard error for Lomax model parameters A, B and B/|A| by
% generating 100 random values from a Lomax distribution with input A and B
% parameters left-truncating/re-zeroing at input theta value and
% right-censoring with input CensorPoint, then optimizing parameters A and B
% individually while holding the other at its known value.
% 
%  Required Inputs:
%      A_opt, B_opt : optimized values for A and B Lomax parameters
%      theta : left-truncation (x-shift) value (used to optimize A_opt and B_opt)
%      CensorPoint : right-censor value 
%      Precision : precision (number of decimal places) for bootstrap sample values
%               randomly generated from lomax distribution with paremeters A_opt and B_opt
% 
%  Optional Inputs:
%      'FitMethod' : 'linear' or 'log' (default = 'linear')
%      'Conditions' : additional parameter conditions (will blow up error if met)
%                (e.g., func_to_opt_pareto(v,X,R,'Condition',{'A>0' 'P<0'}] will NOT allow A>0 or P<0)
%      'Display' : 'on' or 'off' display iterations for fit (default = 'off')
%      'iboot' : number of bootstrap iterations (default = 10000) 
%      'MaxIter' : maximum positive integer number of iterations allowed in fminsearch (default = 100000)
%      'MaxFunEvals' : maximum positive integer number of function evaluations allowed in fminsearch (default = 1000000)
%     
%  Outputs:
%      Ase, Bse, xrse : bootstrapped standard error value for A, B and xr=B/|A| parameters, respectively
%                SE values are standard deviation of bootstrapped A, B and
%                xr values over 
%      A_se, B_se, xr_se :
%  
% 
%  EXAMPLE:
%       
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Parse variable input arguments
okargs={'FitMethod' 'Conditions' 'Display' 'iboot' 'MaxIter' 'MaxFunEvals'};
defaults={'linear' [] 'off' 10000 100000 1000000};
[fitmethod, conditions, display, iboot, maxiter, maxfunevals]=internal.stats.parseArgs(okargs, defaults, varargin{:});

%% Define anonymous functions
col = @(x) reshape(x,[numel(x),1]); %anonymous function to reshape data into single column
DA=@(A,r,c) A(r,c); %select specific elements of array (doesn't require deal) 
trunc=@(x) (x(x>theta)-theta); %anonymous function to left-truncate and re-zero data 


%% Run bootstrap procedure
options = optimset('MaxIter',maxiter, 'MaxFunEvals',maxfunevals,'Display',display); % increase iterations so GP converges

for i=1:iboot

    %Generate 100 values from a Lomax distribution, i.e., Generalized Pareto
    %with theta=0 (input theta value has been truncated, shifted and re-zeroed)
    Xrand = round(random('Generalized Pareto',A_opt, B_opt, 0,[100,1]).*10.^(Precision))./10.^(Precision);
    while any(Xrand<0)
        Xrand(Xrand<0) = round(random('Generalized Pareto',A_opt, B_opt, 0,[size(find(Xrand<0)),1]).*10.^(Precision))./10.^(Precision);
    end
    
    % Calculate empirical distributions for bootstrapped data
    [Xse,~,~,Rse,Pse,~,~] = distributions(trunc(Xrand),'dX',0.001,'CensorPoint',CensorPoint);  
    
    % Optimize lomax parameters for bootstrapped data
    SEcens = Xse>=CensorPoint; 
    [fitresulttest,~,~,~]=createFit(Xse,Rse,'(1+x./lambda).^(-alpha)','fitmethod',fitmethod,'censoring',SEcens, 'Lower',[0 0],'Upper',[inf inf],'Start',[0 0.0001]);
    A_bs_pre0=1./fitresulttest.alpha;
    B_bs_pre= fitresulttest.lambda./fitresulttest.alpha;
    
    [fitresulttestA,~,~,~]=createFit(Xse,Rse,eval(['''(1+A.*x./',num2str(B_bs_pre),').^-(1./A)''']), 'fitmethod',fitmethod,'censoring',SEcens,'upper',[inf],'start',A_bs_pre0,'lower',[-B_bs_pre/max(Xse(~SEcens))],'plot','off');
    A_bs_pre=fitresulttestA.A;
    
    fA = @(v) func_to_opt_pareto(v,Xse,Rse,'censoring',SEcens, 'fitmethod',fitmethod, 'B',B_opt,'theta',0, 'Conditions',conditions, 'AddVar',{'P' Pse});
    A_se(i)=fminsearch(fA,A_bs_pre.*2,options);
    
    fB = @(v) func_to_opt_pareto(v,Xse,Rse, 'censoring',SEcens, 'fitmethod',fitmethod, 'A',A_opt,'theta',0, 'Conditions',conditions, 'AddVar',{'P' Pse});
    B_se(i)=fminsearch(fB,B_bs_pre.*2,options);

    xr_se(i)=B_se(i)./abs(A_se(i));
end    

% Calculate final standard error for A, B and B/|A|
Ase=std(A_se);
Bse=std(B_se);
covAB=DA(cov(col(A_se),col(B_se)),1,2);
xrse=sqrt((B_opt*Ase/A_opt^2)^2 + (Bse/A_opt)^2 - 2*B_opt*covAB/A_opt^3);
