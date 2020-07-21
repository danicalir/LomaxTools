clear
%Load data and set maximum x value
%x1=low_small; x1a=low_small with censored offset; x2=high_large
% load('C:\Users\danicar\Desktop\Analysis\RockDrops\testdata.mat','x1','x2')
% x1=[x1;40];

load('/Users/danicaroth/Google Drive/Academic/Analysis/Scripts/Github/Data_ParticleMotion/Data/rockdrop.mat');
x1=rockdrop{1,1}.data.x;
x2=rockdrop{7,3}.data.x;
D1=0.0168;
x1(x1==0)=D1/2; %make all zero values = D/2
D2=0.0724;
x2(x2==0)=D2/2; %make all zero values = D/2
mindist=0.1;
maxdist=40;
mdist=[mindist maxdist];
xnew=0:0.01:50;

xdata=x1;
D=D1;
xdataa=[xdata;maxdist];

%Calculate PDF and CDF of data and fits
gppl=@(x)[1,1,0;-100,0,0;100,100,D/2];% GP parameter limits for mle [start;lower;upper]
exppl=@(x)[D/2;D/2;max(xdata)];% EXP parameter limits for mle [start;lower;upper]

GP=struct('distname','GeneralizedPareto', 'pdfname','pdf', 'pdf',@(x,k,sigma,theta)gppdf(x,k,sigma,theta), 'cdfname','cdf', 'cdf',@(x,k,sigma,theta)gpcdf(x,k,sigma,theta));
logGP=struct('distname','GeneralizedPareto', 'pdfname','logpdf','pdf',@(x,k,sigma,theta)(-(1+1/k)*log(1+k*(x-theta)./sigma)-log(sigma)), 'cdfname','logsf', 'cdf',@(x,k,sigma,theta)(-(1/k)*log(1+k*(x-theta)/sigma)));
EXP=struct('distname','Exponential', 'pdfname','pdf', 'pdf',@(x,mu)exppdf(x,mu), 'cdfname','cdf', 'cdf',@(x,mu)expcdf(x,mu));
logEXP=struct('distname','Exponential', 'pdfname','logpdf', 'pdf',@(x,mu)((-x/mu)-log(mu)), 'cdfname','logsf', 'cdf',@(x,mu)(-x/mu));




[pd(1)]=fitdistributions(xdata,mdist,GP,gppl(xdata),xnew,'FitMethod','mle');
[pd(2)]=fitdistributions(xdataa,mdist,GP,gppl(xdataa),xnew,'FitMethod','mle');
[pd(3)]=fitdistributions(xdata,mdist,GP,gppl(xdata),xnew,'FitMethod','fitdist','censoring','off');
[pd(4)]=fitdistributions(xdataa,mdist,GP,gppl(xdataa),xnew,'FitMethod','fitdist','censoring','off');
        


% [pd(1)]=fitdistributions(xdata,maxdist,logEXP,exppl(xdata),xnew,'FitMethod','mle');
% [pd(2)]=fitdistributions(xdataa,maxdist,logEXP,exppl(xdataa),xnew,'FitMethod','mle');
% [pd(3)]=fitdistributions(xdata,maxdist,logEXP,exppl(xdata),xnew,'FitMethod','fitdist');%,'censoring','off');
% [pd(4)]=fitdistributions(xdataa,maxdist,logEXP,exppl(xdataa),xnew,'FitMethod','fitdist');%,'censoring','off');

[CDF,X]=ecdf(xdata,'censoring',(xdata>=maxdist));
[CDFa,Xa]=ecdf(xdataa,'censoring',(xdataa>=maxdist));

%Plot everything
% legendtext=['''MLE (censored)'',''MLE (uncensored)'',''Fitdist (censored)'',''Fitdist (uncensored)'''];
legendtext=['''',pd(1).Method,''',''',pd(2).Method,''',''',pd(3).Method,''',''',pd(4).Method,''''];
vtext=@(text) ['xnew,pd(1).',text,',''r'',xnew,pd(2).',text,',''m-.'',xnew,pd(3).',text,',''b'',xnew,pd(4).',text,',''c-.'''];
DC=@(C) deal(C{:});
makevtext=@(varargin) DC(cellfun(@(V) vtext(V),varargin,'uniform',0));
[Ptext,Rtext,PDFtext,CDFtext]=makevtext('P','R','PDF','CDF');

figure(17)
clf
subplot(2,2,1)
hold on
stairs(X,CDF,'k')
stairs(Xa,CDFa,'k-.')
eval(['plot(',CDFtext,')']);
xlabel('x')
xlim([0 max(xdata)]) 
ylabel('Cumulative Probability')
title('Dataset 1')
eval(['legend(''Data'',''Offset Data'',',legendtext,',''Location'',''southeast'')'])

subplot(2,2,2)
[h]=histogram(xdata,100,'Normalization','pdf');
hold on
eval(['plot(',PDFtext,')']);
xlabel('x')
ylabel('Probability Density')
xlim([0 max(xdata)]);
title(pd(1).DistName)

subplot(2,2,3)
eval(['semilogy(',Rtext,')']);
ylabel('Exceedance R')
xlabel('x')
xlim([0 max(xdata)]);

subplot(2,2,4)
eval(['loglog(',Ptext,')']);
ylabel('Disentrainment Rate P')
xlabel('x')
xlim([0 max(xdata)]);
    


