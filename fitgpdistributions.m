function [pd]=fitdistributions(xdata,dist,paramlims,xnew,varargin)

% Fits data to the input distribution with either MLE or fitdist method and
% optional censoring. (Note: fitdist method does not allow censoring for
% generalized Pareto distributions.)
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Required Inputs: 
%   xdata = array of x values
%   dist = structure with fields: 
%       - distname ('GeneralizedPareto', 'Exponential', etc)
%       - pdf and cdf (anonymous functions to call pdf and cdf in mle)
%      e.g., GP=struct('distname','GeneralizedPareto', 'pdf',@(x,k,sigma,theta)gppdf(x,k,sigma,theta), 'cdf',@(x,k,sigma,theta)gpcdf(x,k,sigma,theta));
%            EXP=struct('distname','Exponential', 'pdf',@(x,mu)exppdf(x,mu), 'cdf',@(x,mu) expcdf(x,mu));
%                  
%   paramlims = [starting params;
%              lower bounds;
%              upper bounds]
%      e.g., paramlims_gp=[1,1,0 ; 0,0,0 ; 100,100,min(xdata)];
%   FitMethod = 'mle' or 'fitdist'
% 
% Optional Inputs:
%   censorpoint = [c1,c2] where c1 is lower cutoff and c2 is upper cutoff 
%       - i.e., values of [xdata<c1 | xdata>=c2] will be censored
%       - to skip censoring one or the other, use NaN in place of c1 and/or c2
%   censoring = array of 0 or 1 values indicating elements of xdata to censor
%      e.g., censoring = [xdata<=censorpoint(1) | xdata>=censorpoint(2)];
% 
%  Outputs:
%   pd = data structure with fields:
%         'DistName' : type of distribution
%         'Method' : fitting method
%         'pd' : best fit probabiliity distribution object
%         'X' : x-values at which remaining fields are evaluated
%         'CDF' : model cumulative distribution function values
%         'PDF' : model probability distribution function values
%         'R' : model exceedance probability values
%         'P' : model disentrainment rate values
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%Parse variable input arguments
okargs={'FitMethod' 'CensorPoint' 'Censoring'};
defaults={'fitdist' [] []};
[FitMethod,censorpoint,censoring]=internal.stats.parseArgs(okargs, defaults, varargin{:});

% OPTIONS FOR EACH METHOD AND CALCULATION
vf={'pdf(pd(i).pd,xnew)','cdf(pd(i).pd,xnew)','1-pd(i).CDF','pd(i).PDF./pd(i).R'}; %functions for calculating pdf, cdf, R, P
pdvf={'makedist(dist.distname,param_mle{:})','makedist(dist.distname,param_mle{:})','fitdist(xdata,dist.distname,''censoring'',xdata<=censorpoint(1) | xdata>=censorpoint(2))','fitdist(xdata,dist.distname)'};
valuefunc=[pdvf' [vf;vf;vf;vf]];
valuename={'pd','PDF','CDF','R','P'};
Method={'MLE (censored)','MLE (uncensored)','fitdist (censored)','fitdist (uncensored)'};

% Different approaches for GP and EXP since fitdist won't censor GP
if strcmp(FitMethod,'mle') %Method 1 (MLE)-MAX LIKELIHOOD ESTIMATE FOR GENERALIZED PARETO
    options = statset('MaxIter',10000, 'MaxFunEvals',100000); % increase iterations so GP converges
    if (min(xdata)<=censorpoint(1) || max(xdata)>=censorpoint(2))  && isempty(censoring)==1
        m=1;% MLE Censored data     
        [param_mle,CI_mle]=mle(xdata, dist.pdfname,dist.pdf, dist.cdfname,dist.cdf, 'start',paramlims(1,:), 'lower',paramlims(2,:), 'upper',paramlims(3,:), 'options',options, 'censoring',xdata<censorpoint(1) | xdata>=censorpoint(2)); 
    else
        m=2; % MLE Uncensored data 
        [param_mle,CI_mle]=mle(xdata, dist.pdfname,dist.pdf, dist.cdfname,dist.cdf, 'start',paramlims(1,:), 'lower',paramlims(2,:), 'upper',paramlims(3,:), 'options',options);
    end
    param_mle=num2cell(param_mle);
    CI_mle=num2cell(CI_mle);
    
elseif strcmp(FitMethod,'fitdist') %Method 2 (fitdist)-FITDIST (OR GPFIT SEEMS TO DO THE SAME THING)
    if strcmp(dist.distname,'GeneralizedPareto')
        m=4; %fitdist cannot censor gp
    else
        if (min(xdata)<=censorpoint(1) || max(xdata)>=censorpoint(2)) && isempty(censoring)==1
            m=3;% fitdist Censored data
        else
            m=4; % fitdist Uncensored data
        end
    end
end


for i=1:length(m)
    eval(['pd(i)','=struct(''DistName'',dist.distname,''Method'',Method{m(i)},''pd'',[],''X'',xnew,''CDF'',[],''PDF'',[],''R'',[],''P'',[]);']);
    for j=1:length(valuename)
        eval(['pd(i).',valuename{j},'=',valuefunc{m(i),j},';']);
    end
end