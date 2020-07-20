function [fitresult,gof,fiteq,Rfit] = createFit(X, Y, fitEQ, varargin)

%  Fits input X,Y data to input equation (of any form) by solving
%  for all equation parameters using nonlinear least squares fit method.
%       
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%  Required Inputs:
%      X : independent variable data (1D array)
%      Y : dependent variable data (1D array)
%      fitEQ : model equation to fit for Y (string) (e.g., '(1+x./lambda).^(-alpha)')
%               - Note: string MUST contain independent variable 'x' (lower case), 
%                 and will interpret any other variables as parameters to optimize
%                
%  Optional Inputs:
%      'FitMethod' : 'linear' or 'log' (default = 'linear')
%      'Censoring' : 'X>censorpoint1 | X<censorpoint2' or [elements of X and Y to exclude] (default = [])
%      'Plot' : 'on' or 'off' (default = 'off')
%      'Display' : 'on' or 'off' display iterations for fit (default = 'off')
%      'Lower' : lower limits for [A,B] or [alpha,lambda] (default = [0 0])
%      'Upper' : upper limits for [A,B] or [alpha,lambda] (default = [inf inf])
%      'Start' : starting points for [alpha,lambda] (default = [0 0.0001])
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
okargs={'FitMethod' 'Censoring' 'Plot' 'Display' 'Lower' 'Upper' 'Start'};
defaults={'linear' [] 'off' 'Off' [] [] []};
[fitmethod, censoring, plotstat, display, lower, upper, start]=internal.stats.parseArgs(okargs, defaults, varargin{:});


%% Fit data
% Cut any censored values
Y(censoring)=[];
X(censoring)=[];


% Set up fit type and data depending on whether fit method is log or linear
if strcmpi(fitmethod,'log')   
    [xData, yData] = prepareCurveData( X, log(Y) );
    ft = fittype( ['log(',fitEQ,')'], 'independent', 'x', 'dependent', 'y' );
elseif strcmpi(fitmethod,'linear')
    [xData, yData] = prepareCurveData( X, Y );
    ft = fittype( fitEQ, 'independent', 'x', 'dependent', 'y' );
end


% Set up fit options.
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Display = display;
opts.Lower = lower;
opts.upper = upper;
opts.StartPoint = start;


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


%% Report and display results

% Generate fiteq string for fitting function and results
paramnames = coeffnames(fitresult);
paramvals = coeffvalues(fitresult);
for i=1:length(paramvals)
    eval([paramnames{i},'=paramvals(i);']);
end
eval(['Rfit = @(x) ',fitEQ,';']);
fiteq=['$R = ',fitEQ,'$',newline,'$r^2 = ',num2str(gof.rsquare,'%.2g'),'$'];
paramnames = coeffnames(fitresult);
paramvals = coeffvalues(fitresult);
for i=1:length(paramvals)
    fiteq=[fiteq,newline,'$',paramnames{i},'=',num2str(paramvals(i),'%.2g'),'$'];
end


% Plot fit with data if requested
if strcmpi(plotstat,'on')
    figure(100);
    
    if strcmpi(fitmethod,'log')   
        h = semilogy(X, Y, '.', xData, Rfit(xData), 'r');
    elseif strcmpi(fitmethod,'linear')
        h = plot(X, Y, '.', xData, Rfit(xData), 'r');
    end
    xlabel('X')
    ylabel('R');
    legend( h, 'data', fiteq, 'Location', 'Best','Interpreter','Latex','fontsize',10);
    title([fitmethod,' error minimization'])
end
