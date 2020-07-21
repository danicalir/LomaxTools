function [pd]=fitdistributions(xdata,censorpoint,dist,paramlims,xnew,varargin)

% xdata = array of x values
% censorpoint = [c1,c2] where c1 is lower cutoff and c2 is upper cutoff 
%       - i.e., values of [xdata<=c1 | xdata>=c2] will be censored
%       - to skip censoring, use NaN in place of c1 and/or c2
% dist = structure with fields: 
%       - distname ('GeneralizedPareto', 'Exponential', etc)
%       - pdf and cdf (anonymous functions to call pdf and cdf in mle)
%      e.g., GP=struct('distname','GeneralizedPareto', 'pdf',@(x,k,sigma,theta)gppdf(x,k,sigma,theta), 'cdf',@(x,k,sigma,theta)gpcdf(x,k,sigma,theta));
%            EXP=struct('distname','Exponential', 'pdf',@(x,mu)exppdf(x,mu), 'cdf',@(x,mu) expcdf(x,mu));
%                  
% paramlims = [starting params;
%              lower bounds;
%              upper bounds]
%    e.g., paramlims_gp=[1,1,0 ; 0,0,0 ; 100,100,min(xdata)];
% FitMethod = 'mle' or 'fitdist'
% 

%Parse variable input arguments
okargs={'FitMethod' 'Censoring'};
defaults={[] []};
[FitMethod,censoring]=internal.stats.parseArgs(okargs, defaults, varargin{:});

% OPTIONS FOR EACH METHOD AND CALCULATION
vf={'pdf(pd(i).pd,xnew)','cdf(pd(i).pd,xnew)','1-pd(i).CDF','pd(i).PDF./pd(i).R'}; %functions for calculating pdf, cdf, R, P
pdvf={'makedist(dist.distname,param_mle{:})','makedist(dist.distname,param_mle{:})','fitdist(xdata,dist.distname,''censoring'',xdata<=censorpoint(1) | xdata>=censorpoint(2))','fitdist(xdata,dist.distname)'};
valuefunc=[pdvf' [vf;vf;vf;vf]];
valuename={'pd','PDF','CDF','R','P'};
Method={'MLE (censored)','MLE (uncensored)','fitdist (censored)','fitdist (uncensored)'};

% Different approaches for GP and EXP since fitdist won't censor GP
if length(FitMethod)==3 %Method 1 (MLE)-MAX LIKELIHOOD ESTIMATE FOR GENERALIZED PARETO
    options = statset('MaxIter',6000, 'MaxFunEvals',6000); % increase iterations so GP converges
    if min(xdata)<=censorpoint(1) && max(xdata)>=censorpoint(2)  && isempty(censoring)==1
        m=1;% MLE Censored data     
        [param_mle,CI_mle]=mle(xdata, dist.pdfname,dist.pdf, dist.cdfname,dist.cdf, 'start',paramlims(1,:), 'lower',paramlims(2,:), 'upper',paramlims(3,:), 'options',options, 'censoring',xdata<=censorpoint(1) | xdata>=censorpoint(2)); 
    else
        m=2; % MLE Uncensored data 
        [param_mle,CI_mle]=mle(xdata, dist.pdfname,dist.pdf, dist.cdfname,dist.cdf, 'start',paramlims(1,:), 'lower',paramlims(2,:), 'upper',paramlims(3,:), 'options',options);
    end
    param_mle=num2cell(param_mle);
    CI_mle=num2cell(CI_mle);
    
elseif length(FitMethod)==7 %Method 2 (fitdist)-FITDIST (OR GPFIT SEEMS TO DO THE SAME THING)
    if (min(xdata)<=censorpoint(1) || max(xdata)>=censorpoint(2)) && isempty(censoring)==1
        m=3;% fitdist Censored data
    else
        m=4; % fitdist Uncensored data
    end
end


for i=1:length(m)
    eval(['pd(i)','=struct(''DistName'',dist.distname,''Method'',Method{m(i)},''pd'',[],''X'',xnew,''CDF'',[],''PDF'',[],''R'',[],''P'',[]);']);
    for j=1:length(valuename)
        eval(['pd(i).',valuename{j},'=',valuefunc{m(i),j},';']);
    end
end