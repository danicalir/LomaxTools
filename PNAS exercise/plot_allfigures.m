clear
load('Data/rockdrop.mat');
load('Data/topo.mat');
load('Data/bootstrap.mat');
load('Data/results.mat');%,X,P,R,CensorPoint,A_opt,Ase,B_opt,Bse,xr,xrse,slope,D,site);


%% DEFINE ANONYMOUS FUNCTIONS 
col = @(x) reshape(x,[numel(x),1]); %anonymous function to reshape into columns
DC=@(C) deal(C{:}); %deal elements of cell matrix into separate variables (e.g., outputs of cellfun)
% DA=@(A,r,c) A(r,c); %select specific elements of array (doesn't require deal) 
% DAC=@(A,m,n) deal(A{:}(m,n));
Xfit=@(xmax) 0.001:0.001:xmax;

%% PREP ALL PLOTTING OPTIONS 

iorder=[4 1 5 2 6 3 7]; %order for plotting by slope
[~,q]=sort(iorder); %resorted order for plotting by group (veg/burn)


% Color matrix and rock size delimiter matrix for plotting by site/size
colors=[237 177 32;...
        179 177 32;...
        119 172 48;...
        119 172 48;...
        77 190 238;...
        0 131 205;...
        126 47 142]./255;
plottingsize=(D*1500).^1.3;
linestyle={'-' '--' '-' '--' '-' '--' '-'};
rocklinestyle={'-.' '--' '-'};
eB=1:3;
eV=4:7;
markers=cell(7,3);
markers(eB,:)={'o'};
markers(eV,:)={'^'};
markers=markers(iorder,:);
X_fit=Xfit(50); %continuous array for plotting
eqR=['$R(x) = (1+\frac{A(x-\theta)}{B})^{-1/A}$'];
eqP=['$P(x) = \frac{1}{A(x-\theta)+B}$'];
legtext=['^ $Vegetated$',newline,'o $Burned$',newline,'- 1:1 $line$'];

for s=1:7
    for r=1:3
        dataname{s,r}=[site{s}(1),' (S=',num2str(slopedeg(s,r),2),', D=',num2str(100*D(1,r),2),' cm)'];
    end
end

%% SET PARAMETERS AND RESHAPE VARIABLES FOR PLOTTING

ns=7;nr=3; %numbers of slope and rock variables (for iterations and subplots)
rsvar={X,P,R,CensorPoint,A_opt,Ase,B_opt,Bse,xr,xrse,slope,D,topo.d50}; %cell array of all variables to reorder by slope for plotting
[X,P,R,CensorPoint,A_opt,Ase,B_opt,Bse,xr,xrse,slope,D,d50]=DC(cellfun(@(v)v(iorder,:),rsvar,'Uniform',0)); %reorder all variables
site=site(iorder);
Apos=(A_opt>0); %indices of positive A values
Aneg=(A_opt<0); %indices of negative A values
Pfitopt=@(x,s,r)1./(A_opt(s,r).*(x-theta_opt(s,r))+B_opt(s,r));
Rfitopt=@(x,s,r)(1+(A_opt(s,r).*(x-theta_opt(s,r))./B_opt(s,r))).^(-1./A_opt(s,r));

%% PLOT ROUGHNESS HEIGHT DISTRIBUTIONS (Fig. 2c) 

f2c = figure(1);clf
set(f2c,'Name','Fig. 2c: Cumulative distributions of roughness heights','NumberTitle','off');
hplegtext={};
for i=1:7
    ahp(i)=semilogx(topo.dhp{i},topo.CDF_dhp{i});
    set(ahp(i),'Color',colors(q(i),:),'LineStyle',linestyle{q(i)});
    ylabel('Cumulative fraction less than');
    xlabel('Roughness height d [m]');
    hold on
    
    hplegtext=[hplegtext [site{q(i)}(1:3),', d_{50}=',num2str(topo.d50(i,1),'%.3f')]]; %text for legends to be added later
end
l=legend(hplegtext,'Location','NorthWest');
set(gca,'xlim',[0.0001 3])
set(ahp,'LineWidth',1.25)


%% PLOT ALL MODELED DISENTRAINMENT AGAINST OBSERVED VALUES (Fig. 4)

f4 = figure(3);clf
set(f4,'Name','Fig. 4: Observed vs modeled disentrainment','NumberTitle','off');
for s=1:7
    for r=1:3    
        box on
        scatter((Pfitopt(X{s,r}(X{s,r}<CensorPoint(s,r)),s,r)),(P{s,r}(X{s,r}<CensorPoint(s,r))),plottingsize(s,r),colors(s,:),markers{s,r});
        hold on  
    end
end

set(gca,'yscale','log','xscale','log','xlim',[min(P_fitoptcol) ceil(max(P_fitoptcol)/10)*10],'ylim',[min(Pcol) ceil(max(Pcol)/10)*10])
loglog(get(gca,'xlim'),get(gca,'xlim'),'k'); %one-one line
xlabel('$P_{model}=\frac{1}{Ax+B}$','interpreter','latex','fontsize',20)
ylabel('$P_{observed}$','interpreter','latex','fontsize',20)

% Annotate figure with r-squared
pos = get(gca, 'position');
dim = pos.*[1 1 0.25 0.4]+[.8*pos(3) 0 0 0];
[r,~]=corrcoef(real(log(Pcol)),real(log(P_fitoptcol))); %r-value
rsquare=r(1,2).^2;
annotation('textbox',dim,'String',['$r^2 = ',num2str(rsquare,2),'$'],'LineStyle','none','FaceAlpha',0,'fontsize',14,...
           'FitBoxToText','on','horizontalalignment','left','verticalalignment','bottom','Interpreter','Latex');
       

%% PLOT EXCEEDANCE R VS X WITH OPTIMIZED FITS ON SEMILOG AXES (Fig. 5 insets)

f5R = figure(4);clf
set(f5R,'Name','Fig. 5 (inset): Particle travel distance exceedance with position','NumberTitle','off');
for s=1:7
    for r=1:3     
        % Plot data
        subplot(ns,nr,r+nr*(s-1));
        box on
        hold on
        scatter(X{s,r}(X{s,r}<CensorPoint(s,r)),R{s,r}(X{s,r}<CensorPoint(s,r)),plottingsize(s,r).^.75,colors(s,:),markers{s,r}); %empirical data
        semilogy(Xfit(max(X{s,r}(X{s,r}<CensorPoint(s,r)))),Rfitopt(Xfit(max(X{s,r}(X{s,r}<CensorPoint(s,r)))),s,r),'k','linewidth',1.5); %optimized model line      
        set(gca,'yscale','log','xscale','linear','ylim',[.005 1],'xlim',[0 1.05*max(CensorPoint(s,:))]);
        
        % Title subplots by experiment and label axes
        title(dataname{s,r})
        if s==4 & r==1
            ylabel('exceedance $R$','fontsize',14,'Interpreter','Latex');
        elseif s==7 & r==2
            xlabel('travel distance $x$ (m)','fontsize',14,'Interpreter','Latex')
        end
    end
end


%% PLOT DISENTRAINMENT P VS X WITH OPTIMIZED FITS (Fig. 5 main subplots)

f5P = figure(5);clf
set(f5P,'Name','Fig. 5 (main): Particle disentrainment rates with position','NumberTitle','off');
for s=1:7
    for r=1:3    
        % Plot data
        s5P=subplot(ns,nr,r+nr*(s-1));
        box on
        hold on
        h105a=scatter(X{s,r}(X{s,r}<CensorPoint(s,r)),P{s,r}(X{s,r}<CensorPoint(s,r)),plottingsize(s,r).^.75,colors(s,:),markers{s,r}); %empirical data
        h105c=loglog(Xfit(max(X{s,r}(X{s,r}<CensorPoint(s,r)))),Pfitopt(Xfit(max(X{s,r}(X{s,r}<CensorPoint(s,r)))),s,r),'k','linewidth',1.5); %optimized model line
        set(gca,'yscale','log','xscale','log');
        ylim([0.001 500]);
        xlim([.001 50])
        hold on
        
        % Annotate figures with optimized A, B and B/|A| values +/- errors
        pos = get(s5P, 'position');
        dim = pos.*[1 1 0.5 0.4];
        fitparamsopt{s,r}=['$A = ',num2str(A_opt(s,r),'%.2f'),'\pm',num2str(Ase(s,r),'%.2f'),'$',newline,'$B = ',num2str(B_opt(s,r),'%.2f'),'\pm',num2str(Bse(s,r),'%.2f'),'$',newline,'$\frac{B}{|A|} = ',num2str(abs(xr(s,r)),'%.2f'),'\pm',num2str(xrse(s,r),'%.2f'),'$'];       
        annotation('textbox',dim,'String',fitparamsopt{s,r},'LineStyle','none','FaceAlpha',0,'fontsize',8.5,...
                   'FitBoxToText','on','verticalalignment','bottom','Interpreter','Latex');

        % Title subplots by experiment and label axes
        title(dataname{s,r})     
        if s==4 & r==1
            ylabel('disentrainment rate $P$ [m$^{-1}$]','fontsize',14,'Interpreter','Latex')
        elseif s==7 & r==2
            xlabel('travel distance $x$ (m)','fontsize',14,'Interpreter','Latex')
        end
    end
end


%% PLOT LOMAX PARAMETERS AGAINST PHYSICAL VARIABLES (Fig. 6)

f6 = figure(6);clf
set(f6,'Name','Fig. 6: Lomax parameters','NumberTitle','off');

% Set cell arrays with arguments to use in for loop
LomaxParamname={{'$\frac{B}{|A|}$'} {'$|A|$'} {'$B$'}}; %parameter strings
LomaxParam={abs(B_opt./A_opt) abs(A_opt)  B_opt}; %parameter arrays
LomaxErr={xrse Ase Bse}; %error arrays
Lomaxfunc={ @(a)(a) @(a)abs(a) @(a)(a)}; %anonymous functions for calculating
xvalname={'$SD/d_{50}$' '$SD$'}; %x-axis strings
xval={slope.*D./d50 slope.*D}; %x-axis arrays
axisscale={'''xscale'',''log'',''yscale'',''log'''}; 
corrcoefarg={{'real(log(col(xval{j}(Apos)))),real(log(col(LomaxParam{i}(Apos))))'};...
             {'real(log(col(xval{j}(Aneg)))),real(log(col(LomaxParam{i}(Aneg))))'};...
             {'real(log(col(xval{j}))),real(log(col(LomaxParam{i})))'}};  
corrcoefname={'_+' '_-' ''};
corrcoefselect={{1:2 1:2} {1:2 1:2} {1:2 3}};
corrcoefpos={[0 .8]};
xlimit=cellfun(@(a) [.9*min(a,[],'all') 1.1*max(a,[],'all')],xval,'uniformoutput',0);
negmarkershape=markers;
for i=1:length(LomaxParam) 
    for j=1:length(xval)
        for s=1:7
            for r=1:3
                subplot(length(LomaxParam),length(xval),j+length(xval)*(i-1))
                box on
                if A_opt(s,r)<0
                    scatter(col(xval{j}(s,r)),col(Lomaxfunc{i}(LomaxParam{i}(s,r))),plottingsize(s,r),colors(s,:),negmarkershape{s,r},'filled');
                end 
                hold on
                if abs(LomaxParam{i}(s,r))<abs(LomaxErr{i}(s,r))
                    errorbar(col(xval{j}(s,r)),col(Lomaxfunc{i}(LomaxParam{i}(s,r))),col(Lomaxfunc{i}(LomaxParam{i}(s,r)))-.01,...
                        NaN,'color','k','LineStyle','-.','marker','none','LineWidth',.5)
                    errorbar(col(xval{j}(s,r)),col(Lomaxfunc{i}(LomaxParam{i}(s,r))),NaN,col(Lomaxfunc{i}(LomaxErr{i}(s,r))),...
                        'color',colors(s,:),'marker','none','LineWidth',.5)
                else 
                    errorbar(col(xval{j}(s,r)),col(Lomaxfunc{i}(LomaxParam{i}(s,r))),col(Lomaxfunc{i}(LomaxErr{i}(s,r))),'color',colors(s,:),'marker','none','LineWidth',.5)
                end
                scatter(col(xval{j}(s,r)),col(Lomaxfunc{i}(LomaxParam{i}(s,r))),plottingsize(s,r),colors(s,:),'marker',markers{s,r},'LineWidth',1.5);
                if i==3
                    xlabel(xvalname{j},'Interpreter','Latex','fontsize',14)
                end
            end 
        end
        xlim(xlimit{j});
        eval(['set(gca,',axisscale{:},')']); 
        
        if j==1
            ylabel(LomaxParamname{i},'Interpreter','Latex','fontsize',14)
        end
        
        if ismember(j,corrcoefselect{i}{1})
            rptext=[];
            for k=1:length(corrcoefselect{i}{2})
                rval=eval(['corrcoef(',corrcoefarg{corrcoefselect{i}{2}(k)}{:},');']);
                rptext=[rptext,newline,'$r',corrcoefname{corrcoefselect{i}{2}(k)},' = ',num2str(rval(1,2),2),'$'];
                        
            end
            pos = get(gca, 'position');
            dim = [pos(1)+.7*pos(3) pos(2)+.1*pos(4) 0.1 0.1];
            dim = [pos(1)+corrcoefpos{1}(1)*pos(3) pos(2)+corrcoefpos{1}(2)*pos(4) 0.5 0.4];
            annotation('textbox',dim,'String',rptext,'LineStyle','none','FaceAlpha',0,'fontsize',11,...
                       'FitBoxToText','on','verticalalignment','bottom','Interpreter','Latex');
        end
    end
end

pause(2)
xlim(xlimit{j})