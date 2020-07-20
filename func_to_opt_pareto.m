function SSE = func_to_opt_pareto(v,X,R,varargin)

% Calcuates the sum of squared errors (SSE) between input X and R data arrays
% and a generalized Pareto (or Pareto II) model distribution (form below) with
% parameters A, B and theta.
% 
%   Generalized Pareto distribution:
% 
%        R = [1 + A(X-theta)/B]^(-1/A),  for A =/= 0
%            exp[-(X-theta)/B]           for A = 0.
% 
% Note: func_to_opt can be used as part of an joint optimization algorithm for
%       A, B and theta by minimizing SSE with fminsearch (see example below).
%        
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%  Required Inputs:
%      X: 1D independent variable array 
%      R: 1D exceedance probability array 
%               (R = 1 - CDF, where CDF is the empirical cumulative distribution)
%      v: vector containing [A B theta] values 
%               For use with the fminsearch optimization algorithm, v is an
%               anonymous function input used optimize A, B and theta values 
%               (leave as 'v'; see example below)  
% 
%  Optional Inputs:
%      'FitMethod' : 'linear' or 'log' values of errors to minimize (default = 'linear')
%      'Censoring' : 'X>censorpoint1 | X<censorpoint2' or [elements of X to exclude] (default = [])
%      'Conditions' : additional parameter conditions (will blow up error if met)
%                (e.g., func_to_opt(v,X,R,'Conditions',{'A>0' 'P<0'}] will NOT allow A>0 or P>0)
%      'AddVars' : cell array of additional variables needed to evaluate added Conditions
%                - must be in format {'variable1name' variable1; 'variable2name' variable2;...}
%                (e.g., to include P from above example condition, use: 'AddVars',{'P' P})
%      'A' : set value for A
%      'B' : set value for B
%      'theta' : set value for theta
% 
% 
% EXAMPLE: Use func_to_opt_pareto as part of an optimization algorithm for
%          A, B and theta by minimizing SSE with fminsearch: 
% 
%    f = @(v) func_to_opt(v,X,R,'censoring',cens, 'fitmethod','linear', 'theta',0, 'Conditions',conditions);
%    vopt = fminsearch(f,[A,B,theta]); 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Parse variable input arguments
okargs={'FitMethod' 'Censoring' 'Conditions' 'AddVars' 'A' 'B' 'theta'};
defaults={'linear' [] [] [] [] [] [] []};
[fitmethod, censoring, addconditions, addvars, A, B, theta]=internal.stats.parseArgs(okargs, defaults, varargin{:});

% Censor values if any entered
R(censoring)=[];
X(censoring)=[];


%% Define model, parameters and outputs

% Define model parameters to optimize
ABtheta=[isempty(A) isempty(B) isempty(theta)];
ABthetai=find(ABtheta==1);
ABthetalabels = {'A' 'B' 'theta'};
for i=1:sum(ABtheta)
    eval([ABthetalabels{ABthetai(i)},'=v(i);']);
end

% Define model
if A~=0
    Rhat = (A.*(X-theta)./B+1).^(-1./A);
    Xhat = (((R.^(-A))-1).*B./A)+theta;

elseif A==0
    Rhat = exp(-X./B);
    Xhat = -B.*log(R);
end  

% Define residuals (log or linear model) to use for minimization
if strcmpi(fitmethod,'log')   
    errsR = log(Rhat)-log(R);
    errsX = log(Xhat)-log(X);
    
elseif strcmpi(fitmethod,'linear') 
    errsR = (Rhat)-(R);
    errsX = (Xhat)-(X);
end

% Define sum of squared errors SSE
    SSE = sum(errsR.^2); 

    
    
%% Set forbidden parameter conditions by blowing up error 
if (B-A*theta) < 0
    SSE = 1e20;
end

if A<0 && max(X)>(A*theta-B)/A
    SSE = 1e20;
end

if A>0 && 0<(A*theta-B)/A
    SSE = 1e20;
end


% Add any additional variables (needed for parameter condition requirements)
for i=1:size(addvars,1)
    eval([addvars{i,1},'=addvars{i,2};']);
end

% Apply any additional parameter condition requirements
for i=1:length(addconditions)
    if eval(addconditions{i})
        SSE = 1e20;
    end
end