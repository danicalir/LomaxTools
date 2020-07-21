% Compares Lomax optimization using lomaxopt versus built-in MLE and fitdist (which produces
% same output as gpfit) functions
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear

load('Data/rockdrop.mat');
for s=1:7
    for r=1:3
theta=0.1; %left truncation point (m)
trunc=@(data) (data.x(data.x>theta)-theta); %anonymous function for left-truncating and re-zeroing data
maxdist=rockdrop{s,r}.distributions.CensorPoint;
mdist=[0 maxdist];
xnew=0:0.001:maxdist;

xdata=trunc(rockdrop{s,r}.data);

%Calculate PDF and CDF of data and fits
gppl=[1,1,0;-100,0,0;100,100,0];% GP parameter limits for mle [start;lower;upper]

GP=struct('distname','GeneralizedPareto', 'pdfname','pdf', 'pdf',@(x,k,sigma,theta)gppdf(x,k,sigma,theta), 'cdfname','cdf', 'cdf',@(x,k,sigma,theta)gpcdf(x,k,sigma,theta));
logGP=struct('distname','GeneralizedPareto', 'pdfname','logpdf','pdf',@(x,k,sigma,theta)(-(1+1/k)*log(1+k*(x-theta)./sigma)-log(sigma)), 'cdfname','logsf', 'cdf',@(x,k,sigma,theta)(-(1/k)*log(1+k*(x-theta)/sigma)));

[pd(1)]=fitgpdistributions(xdata,GP,gppl,xnew,'FitMethod','mle','CensorPoint',mdist);
[pd(2)]=fitgpdistributions(xdata,GP,gppl,xnew,'FitMethod','fitdist','CensorPoint',mdist);

cens=rockdrop{s,r}.distributions.X>=rockdrop{s,r}.distributions.CensorPoint; % Elements of X to censor from Lomax fitting

[A_opt,B_opt] = lomaxopt(rockdrop{s,r}.distributions.X,rockdrop{s,r}.distributions.R, 'censoring',cens, 'fitmethod','log' , 'conditions',{'any(Rhat<1e-6)'  '1./(B-A.*theta)<0'} , 'addvars',{'P' rockdrop{s,r}.distributions.P});

if any(xdata>=maxdist)
    censmethod='(censored)';
elseif ~any(xdata>=maxdist)
    censmethod='(uncensored)';
end

%Plot everything
legendtext={['Data ',censmethod,''] ['',pd(1).Method,''] ['',pd(2).Method,''] ['lomaxopt ',censmethod,'']};
vtext=@(text) ['xnew,pd(1).',text,',''b'',xnew,pd(2).',text,',''m-.'',xnew,',text,'fitopt(xnew),''g--'''];%,xnew,pd(3).',text,',''b'',xnew,pd(4).',text,',''c-.'''];
DC=@(C) deal(C{:});
makevtext=@(varargin) DC(cellfun(@(V) vtext(V),varargin,'uniform',0));
[Ptext,Rtext,PDFtext,CDFtext]=makevtext('P','R','PDF','CDF');
Pfitopt=@(x)1./(A_opt.*x+B_opt);
Rfitopt=@(x)(1+(A_opt.*x./B_opt)).^(-1./A_opt);
PDFfitopt=@(x)(1/B_opt).*(1+(A_opt.*x./B_opt)).^(-((1./A_opt)+1));
CDFfitopt=@(x)(1-Rfitopt(x));
xrange=[0 max(xdata).*1.1];

figure(17)
clf
subplot(2,2,1)
hold on
stairs(rockdrop{s,r}.distributions.X,rockdrop{s,r}.distributions.CDF,'k')
eval(['plot(',CDFtext,')']);
xlabel('x')
xlim(xrange) 
ylabel('Cumulative Probability')
title(rockdrop{s,r}.data.SiteName)
legend(legendtext,'Location','southeast');

subplot(2,2,2)
[h]=histogram(xdata,100,'Normalization','pdf');
hold on
eval(['plot(',PDFtext,')']);
xlabel('x')
ylabel('Probability Density')
xlim(xrange) 

subplot(2,2,3)
eval(['semilogy(rockdrop{s,r}.distributions.X,rockdrop{s,r}.distributions.R,''k.'',',Rtext,')']);
ylabel('Exceedance R')
xlabel('x')
xlim(xrange) 

subplot(2,2,4)
eval(['loglog(rockdrop{s,r}.distributions.X,rockdrop{s,r}.distributions.P,''k.'',',Ptext,')']);
ylabel('Disentrainment Rate P')
xlabel('x')
xlim(xrange) 
   
pause

    end
end

