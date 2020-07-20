function [X,f,F,R,P,Xemp,Femp,CensorPoint] = distributions(data,varargin)

% GENERATES EMPIRICAL DISTRIBUTIONS FOR EMPIRICAL DATA VECTOR WITH CENSORING
% 
% INPUT: - data: data structure with fields:
%               [.x (data vector)] OR [.CDF and .X (e.g., empirical CDF from ecdf)]  
%               [.Censorpoint (or manual censorpoint may be input)]
%               [.IntermediateAxis (or manual starting point may be input)]
%          OR 1D vector of sample values (e.g., [nx1] vector)
%
% OPTIONAL INPUT:
%         - dX: dX used for calculating empirical pdf 
%               - default dX = 10^(-p), where p is the highest precision of any x measurement
%                 (e.g., sample array x = [1, 2.17, 3.8] would lead to dX = 0.01)
%         - CensorPoint: scalar right censor value with same units as x or X (default no censorpoint)
%               - may also be included as .CensorPoint field in data (but
%                 manually input 'CensorPoint' argument will override)
%         - offset: 'on' or 'off' (default on)
%               - 'on' assumes 100% of probability is NOT contained within the sampled
%                 range [x<CensorPoint] even if 100% of sampled values are <CensorPoint,
%                 ie, uncensored CDF < 1 at maximum sample value 
%                 (equivalent to calculating distributions using 1/(N+1) or 1/(N+w) normalization)
%               - 'off' assumes samples adequately characterize 100% of probability within the system,
%                  ie, uncensored CDF = 1 at maximum sample value 
%                 (equivalent to using 1/N normalization)
%         - offsetweight: numerical value w used to weight normalization offset (default 1)
%               - for data.x or x-vector (data) input, interpreted as integer number of points: 
%                   x will be weighted with w additional points at x > CensorPoint
%                   (equivalent to 1/(N+w) normalization)
%               - for data.CDF (CDF) input, interpreted as percent probability reduction for renormalization:
%                   CDF will be renormalized by multiplying by (1-w/100),
%                   ie, uncensored CDF = 1-w/100 at maximum sample value
%                   (equivalent to (1/N)*(1-w/100) normalization)
%
% OUTPUT: - X: censored data vector
%         - f: empirical probability density distribution
%         - F: empirical cumulative probability distribution
%         - R: empirical survival function
%         - P: empirical disentrainment rate
%         - Xemp: full empirical data vector
%         - Femp: full empirical CDF
%         - CensorPoint: value used to right-censor distributions
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% Written by Danica Roth, University of Oregon, October 2018.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%Parse variable input arguments
okargs={'dX' 'CensorPoint' 'offset' 'offsetweight'};
defaults={[] NaN 'on' 1};
[dX,CensorPoint,offset,offsetweight]=internal.stats.parseArgs(okargs, defaults, varargin{:});

msg = 'Data must be input as either 1D sample vector or data structure with fields ''x'' OR ''CDF'' and ''X''';
%If data is input as sample array, reformat as data structure
if ~isstruct(data)
    if isvector(data)
        data=struct('x',data);
    else
        error(msg);
    end
end

% If input is data array, reshape it for additional processing
if isfield(data,'x')
    x=reshape(data.x,[length(data.x),1]);
    
    % If offset is turned on, add number of extra points indicated by offsetweight past Censor value
    if strcmpi(offset,'on')
        x=[x;repmat(CensorPoint+10,offsetweight,1)];
    end

    %Empirical CDF from data
    [Femp,Xemp]=ecdf(x); %must include any censored points or interpolation can drop the endpoint

% If input is CDF, reshape it
elseif isfield(data,{'CDF','X'})  
    Femp=reshape(data.CDF,[length(data.CDF),1]);
    Xemp=reshape(data.X,[length(data.X),1]);
    
    % If offset is turned on, renormalize
    if strcmpi(offset,'on')
        Femp=Femp.*(1-offsetweight/100);
    end
else
    error(msg);
end

% Set CensorPoint 
if isfield(data,'CensorPoint') & isnan(CensorPoint) %if data includes .CensorPoint field and no value is set on input, use field value
    CensorPoint=data.CensorPoint;
elseif ~isfield(data,'CensorPoint') & isnan(CensorPoint)  %if data does not includes .CensorPoint field and no value is set on input, use 110% max data value (ie, do not censor)
    CensorPoint=1.1*max(Xemp);
end


%Spacing for interpolation and empirical pdf calculation
%If no dX is specified, use precision of x-values to select dX
if isempty(dX)
    dX=min(arrayfun(@(x) getprecision(x),Xemp(Xemp>0)));
end

%Evenly spaced X vector for interpolation
X=[0:dX:(max(Xemp(Xemp<CensorPoint))+dX)]';

%Add zero for interpolation at edges in case probability starts at Xemp>0 (gind will remove duplicates)
Femp=[0;Femp];
Xemp=[0;Xemp];
gind=[(eldiff(Xemp)~=0);true]; %ecdf creates a duplicate point by assuming F=0 at first X value (remove this)
Femp=Femp(gind);
Xemp=Xemp(gind);

%Calculate evenly interpolated probability distributions
F=interp1(Xemp,Femp,X,'linear',max(Femp)); %interpolate CDF for equal dX only up to max(X)<CensorPoint; extrapolate with final value to avoid NaNs
X(isnan(F))=[]; %get rid of any edge values for which F is NaN 
F(isnan(F))=[]; 
f=diff([0;F])./dX; %empirical probability density distribution

%Keep final arrays of only observed Xemp values < CensorPoint
cen=find(Xemp<CensorPoint,1,'last');
F=interp1(X,F,Xemp(2:cen));
f=interp1(X,f,Xemp(2:cen));
X=Xemp(2:cen); %rename final observed x-array for function output
R=1-F; %empirical survival function
P=f./R; %empirical disentrainment function    
end


function p = getprecision(x)
f = 14-9*isa(x,'single'); % 14 (double) or 5 (single) decimal places 
s = sprintf('%.*e',f,x); % full value as string
v = [f+2:-1:3,1];
s(v) = '0'+diff([0,cumsum(s(v)~='0')>0]);
p = str2double(s);
end