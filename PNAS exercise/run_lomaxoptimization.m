% Run lomax optimization for all sites and particle sizes.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear
load('Data/rockdrop.mat');
load('Data/topo.mat');

%% DEFINE ANONYMOUS FUNCTIONS
% col = @(x) reshape(x,[numel(x),1]); %anonymous function to reshape into columns
% DC=@(C) deal(C{:}); %deal elements of cell matrix into separate variables (e.g., outputs of cellfun)
% DA=@(A,r,c) A(r,c); %select specific elements of array (doesn't require deal) 
% DAC=@(A,mi) deal(A{:}(mi));
% Xfit=@(xmax) 0.001:0.001:xmax;



%% SET PARAMETERS AND CALCULATE VARIABLES 
%numbers of slope and rock variables (for iterations and subplots)
ns=7;
nr=3;
% Pull parameters from 'rockdrop' structure
% rockdrop=rockdrop(iorder,:);
site=cellfun(@(v)v(1:3),topo.site,'uniform',0); %dataset site identifiers (first letter of each site filename: H=HPB, topo=Noble)
slopedeg=repmat(topo.slopedeg',1,3);%hillslope angles in degrees
sloperad=deg2rad(slopedeg); %hillslope angles converted to radians
slope=repmat(topo.slope',1,3); %hillslope gradient
D=repmat(cellfun(@(c)c.data.IntermediateAxis,rockdrop(1,:)),7,1); %rock diameters [m]
Fx=cellfun(@(v)v.distributions.CDF,rockdrop,'Uniform',0); %CDF values associated with X array
P=cellfun(@(v)v.distributions.P,rockdrop,'Uniform',0); %Disentrainment rates associated with X array
R=cellfun(@(v)v.distributions.R,rockdrop,'Uniform',0); %Exceedance associated with X array
X=cellfun(@(v)v.distributions.X,rockdrop,'Uniform',0); %Unique particle travel distances 
theta=0.1; %left truncation point (m) (shift in x)
CensorPoint=cell2mat(cellfun(@(v)v.distributions.CensorPoint,rockdrop,'Uniform',0));

fitmethod='log';
conditions={'any(Rhat<1e-6)'  '1./(B-A.*theta)<0'};

% Make columns
Xcol=[];
Rcol=[];
Pcol=[];
Scol=[];
Dcol=[];
CensorPointcol=[];
for s=1:size(X,1)
    for r=1:size(X,2)
        N(s,r)=numel(X{s,r});
        Xcol=[Xcol;X{s,r}];
        Rcol=[Rcol;R{s,r}];
        Pcol=[Pcol;P{s,r}];
        Scol=[Scol;repmat(slope(s,r),N(s,r),1)];
        Dcol=[Dcol;repmat(D(s,r),N(s,r),1)];
        CensorPointcol=[CensorPointcol;repmat(CensorPoint(s,r),N(s,r),1)];
    end
end


%% OPTIMIZE LOMAX FITS (A,B) FOR P AND R 
P_fitprecol=[];
R_fitprecol=[];
Acol=[];
Bcol=[];
thetacol=[];
P_fitoptcol=[];
R_fitoptcol=[];
for s=1:7 %slope identifier    
    for r=1:3 %rock size identifier
        cens=X{s,r}>=CensorPoint(s,r); % Elements of X to censor from Lomax fitting
        [A_opt(s,r),B_opt(s,r)] = lomaxopt(X{s,r},R{s,r}, 'censoring',cens, 'fitmethod',fitmethod, 'conditions',conditions, 'addvars',{'P' P{s,r}});
        
        theta_opt(s,r)=0;
        xr(s,r)=B_opt(s,r)./abs(A_opt(s,r));
        
        Rfitopt=@(x,s,r)(1+(A_opt(s,r).*(x-theta_opt(s,r))./B_opt(s,r))).^(-1./A_opt(s,r));
        Pfitopt=@(x,s,r) 1./(A_opt(s,r).*(x-theta_opt(s,r))+B_opt(s,r));
        
        % Compile columns
        Acol=[Acol;repmat(A_opt(s,r),N(s,r),1)];
        Bcol=[Bcol;repmat(B_opt(s,r),N(s,r),1)];
        thetacol=[thetacol;repmat(theta_opt(s,r),N(s,r),1)];
        
        P_fitopt{s,r}=Pfitopt(X{s,r},s,r);
        P_fitoptcol=[P_fitoptcol;P_fitopt{s,r}];

        R_fitopt{s,r}=Rfitopt(X{s,r},s,r);
        R_fitoptcol=[R_fitoptcol;R_fitopt{s,r}];
     end
end 


%% Cut model exceedance values if R<0.001 (estimated error tolerance threshold) or R>1 (probability limit).
Routind=find(R_fitoptcol<0.001 | R_fitoptcol>1);
Xoutind=find(Xcol<=theta | Xcol>=CensorPointcol);
outind=[Routind;Xoutind];
R_fitout=R_fitoptcol(outind);
P_fitout=P_fitoptcol(outind);
Xout=Xcol(outind);
Rout=Rcol(outind);
Pout=Pcol(outind);
Xcol(outind)=[];
R_fitoptcol(outind)=[];
P_fitoptcol(outind)=[];
Rcol(outind)=[];
Pcol(outind)=[];
Acol(outind)=[];
Bcol(outind)=[];
thetacol(outind)=[];
Dcol(outind)=[];
Scol(outind)=[];


save('Data/results.mat');

disp('Results of lomax optimization saved to Data/results.mat');




