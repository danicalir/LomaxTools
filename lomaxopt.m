function [A_opt,B_opt] = lomaxopt(X,R,varargin)

%  Optimizes parameters A and B for input X and R data fit to a Lomax
%  distribution (Generalized Pareto or Pareto II with location parameter
%  theta set to 0) defined by:
%       
%  Lomax distribution:
% 
%        R = [1 + AX/B]^(-1/A),  for A =/= 0
%            exp(-X/B)           for A = 0.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   
%  Required Inputs:
%      X: X data 
%      R: R data, R = [1+A(x-theta)/B]^-1/A or [1+(x-theta)./lambda].^(-alpha)
% 
%  Optional Inputs:
%      'FitMethod' : 'linear' or 'log' (default = 'log')
%      'Censoring' : 'X>censorpoint1 | X<censorpoint2' or [elements of X and Y to exclude] (default = [])
%      'Conditions' : additional parameter conditions (will blow up error if met)
%                (e.g., func_to_opt_pareto(v,X,R,'Conditions',{'A>0' 'P<0'}] will NOT allow A>0 or P<0)
%      'AddVars' : cell array of additional variables needed to evaluate added Conditions
%                - must be in format {'variable1name' variable1; 'variable2name' variable2;...}
%                (e.g., to include P from above example condition, use: 'AddVars',{'P' P})
%      'MaxIter' : maximum positive integer number of iterations allowed in fminsearch (default = 100000)
%      'MaxFunEvals' : maximum positive integer number of function evaluations allowed in fminsearch (default = 1000000)
%      'Display' : 'on' or 'off' display iterations for fit (default = 'off')
% 
%  Outputs:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%      fiteq : string reporting fit equation, r-squared, parameter names and values 
% 
% 
% EXAMPLE:
% [fitresult,gof,fitparams,Rfit] = createFit(X,R,'fitmethod',fitmethod,'censoring',cens,'upper',[inf],'start',[0 .01],...
%                                            'fitEQ','(1+A.*x./B).^-(1./A)','lower',[],'plot','off');  %
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



%% Parse variable input arguments
okargs={'FitMethod' 'Censoring' 'Conditions' 'AddVars' 'MaxIter' 'MaxFunEvals' 'Display'};
defaults={'linear' [] [] [] 100000 1000000 'off'};
[fitmethod, cens, conditions, addvars, maxiter, maxfunevals, display]=internal.stats.parseArgs(okargs, defaults, varargin{:});

% Read in any extra variables added
if ~isempty(addvars)
    for i=1:size(addvars,1)
        eval([addvars{i,1},'=addvars{i,2};']);
    end
end

%% Fit preliminary Lomax parameters
% (used to set reasonable starting values for dual optimization algorithm below)

% Use equivalent Pareto II distribution form with A=1/alpha and B=lambda/alpha:
% R = [1+(x-theta)./lambda].^(-alpha) to avoid 

[fitresultR,~,~,~]=createFit(X,R,['(1+x./lambda).^(-alpha)'],'fitmethod',fitmethod,'censoring',cens,...
                                           'lower',[0 0],'upper',[inf inf],'start',[0 0.0001]);

A_pre0=1./fitresultR.alpha; %preliminary A value used as initial starting value for optimizing A
B_pre=fitresultR.lambda./fitresultR.alpha; %preliminary B value, then set and fit for A

%A fit to B set from R-fit (setting B allows constraints on A so negative A possible)
[fitresultRA, ~,~,~] = createFit(X,R,eval(['''(1+A.*x./',num2str(B_pre),').^-(1./A)''']),'fitmethod',fitmethod,...
                                           'censoring',cens,'upper',[inf],'start',A_pre0,'lower',[-B_pre/max(X(~cens))],'plot','off');
A_pre = fitresultRA.A;


%% Use fminsearch for joint optimization of A and B 
% (with starting values from above fitting procedure)

f = @(v) func_to_opt_pareto(v,X,R,'censoring',cens, 'fitmethod',fitmethod, 'theta',0, 'Conditions',conditions, 'AddVars',{'P' P});
options = optimset('MaxIter',maxiter, 'MaxFunEvals',maxfunevals, 'Display',display); % increase iterations so GP converges

vopt = fminsearch(f,[A_pre.*2 B_pre.*2],options); %use previous rough fits as starting values

A_opt=vopt(1);
B_opt=vopt(2);
