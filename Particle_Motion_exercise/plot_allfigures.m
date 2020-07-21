

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

fc=[     0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];%colors for histograms  
yt=[300;...%transect y(element) number (manually selected for each site)
    150;...
    120;...
    94;...
    220;...
    155;...
    304];
yt=yt(iorder);
isite='BBBNNNN';
hnum=[10 30 30]; %number of histogram bins for different particle sizes

%% SET PARAMETERS AND RESHAPE VARIABLES FOR PLOTTING

ns=7;nr=3; %numbers of slope and rock variables (for iterations and subplots)
rsvar={X,P,R,CensorPoint,A_opt,Ase,B_opt,Bse,xr,xrse,slope,slopedeg,D,topo.d50}; %cell array of all variables to reorder by slope for plotting
[X,P,R,CensorPoint,A_opt,Ase,B_opt,Bse,xr,xrse,slope,slopedeg,D,d50]=DC(cellfun(@(v)v(iorder,:),rsvar,'Uniform',0)); %reorder all variables
site=site(iorder);
Apos=(A_opt>0); %indices of positive A values
Aneg=(A_opt<0); %indices of negative A values
Pfitopt=@(x,s,r)1./(A_opt(s,r).*(x-theta_opt(s,r))+B_opt(s,r));
Rfitopt=@(x,s,r)(1+(A_opt(s,r).*(x-theta_opt(s,r))./B_opt(s,r))).^(-1./A_opt(s,r));
theta=0.1; %theta (left truncation) at all sites = 0.1 m

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


%% PLOT TOPOGRAPHIC PROFILES AND TRAVEL DISTANCE HISTOGRAMS (Fig. 2a, 2b and S2)

for i = 1:7
    eval(['f2',site{i},' = figure(2',num2str(i),');']); clf
    eval(['set(f2',site{i},',''Name'',''Figs. 2/S2: Topography and travel distance histograms'',''NumberTitle'',''off'');']);

    %Pick rasters
    raster={topo.raster{iorder(i)} topo.raster_lp{iorder(i)} topo.raster_hp{iorder(i)}};
    nx=size(raster{1},2); 
    ny=size(raster{1},1);
    y = [1:ny]./100;   
    x = [1:nx]./100;
    z = raster{1}(yt(i),:);
    zs = z-(slope(i).*x); %re-trend elevations in horizontal reference frame
    x = x./cosd(slopedeg(i))-theta; %convert x to along-ground reference frame (matching travel distances) and shift by truncation distance (theta)
    z0 = interp1(x(~isnan(zs)),zs(~isnan(zs)),0,'linear','extrap'); %extrapolate linearly to find z(x=0) to rezero z (for sites with z(x=0)=NaN)
    zs = zs-z0; %rezero z to x=0 (theta-shifted starting line) elevation
    cp(i)=max(CensorPoint(i,:)-theta); %make i-specific cp and shift by theta for easy reference later
       
    %Plot the topographic transects 
    s1 = subplot(2,1,1);
    yyaxis left
    plot(x,zs,'k','linewidth',2); hold on
    ylabel('elevation below starting point [m]')
    
    %Set lower z-limit for plotting
    if i==1
        zlimval = [-11 0]; %manually set first site because raster does not extend to max travel distance
    else
        zlimval = [interp1(x(~isnan(zs)),zs(~isnan(zs)),1.05*min(max(x),cp(i)),'linear','extrap') 0];
    end
    
    %Plot the filtered topography
    yyaxis right
    plot(x,raster{2}(yt(i),:),'-s','linewidth',5,'Color',0.85.*[1 1 1]); %plot low-pass (long wavelength) topographic transects (in gray for visibility)
    plot(x,raster{3}(yt(i),:),'linewidth',2,'Color',colors(i,:)); %plot high-pass (short wavelength) topographic transects)
    ylabel('filtered topography [m]')    
    set(s1.YAxis(1),'Color','k','Limits',zlimval)
    set(s1.YAxis(2),'Color',colors(i,:),'Limits',[-.75 .75])
    set(s1.XAxis,'Limits',[0 1.05*cp(i)]);
    title(site{i}(1:3));
    
    %Plot close-up inset
    axes('Position',[.64 .84 .25 .07])
    box on
    plot(x,raster{3}(yt(i),:),'linewidth',2,'Color',colors(i,:));
    ylabel('d [m]')
    xlabel('x [m]')
    ylim([-.3 .3])
    xlim([0 2.3])
    ax2=gca;
    ax2.YAxis(1).Color = colors(i,:);

    %Plot histograms of travel distances
    subplot(2,1,2)
    for r=[3 2 1]
        xdata{r} = rockdrop{iorder(i),r}.data.x(rockdrop{iorder(i),r}.data.x>0.1)-0.1;
        xdata{r}(xdata{r}>=cp(i))=cp(i);
        histogram(xdata{r},hnum(r),'FaceColor',fc(r,:)); hold on
    end
    ylabel('Number of particles disentrained')
    xlabel('along-ground travel distance x [m]')
    xlim([0 1.05*cp(i)]) %match axes with final figure 2 (x vs R)
  
    if i==1
        legend('large','medium','small')
    end
end


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
corrcoefarg={{'log(col(xval{j}(Apos))), log(abs(col(LomaxParam{i}(Apos))))'};...
             {'log(col(xval{j}(Aneg))),  log(abs(col(LomaxParam{i}(Aneg))))'};...
             {'log(col(xval{j})), log(abs(col(LomaxParam{i})))'}};  
corrcoefname={'_+' '_-' ''};
corrcoefselect={{1:2 1:2} {1:2 1:2} {1:2 3}};
corrcoefpos={[0 .8]};
xlimit=cellfun(@(a) [.9*min(a,[],'all') 1.1*max(a,[],'all')],xval,'uniformoutput',0);
ylimit={[1e-2 1e3] [1e-3 1e1] [1e-2 1e1]};

for i=1:length(LomaxParam) 
    for j=1:length(xval)
        for s=1:7
            for r=1:3
                subplot(length(LomaxParam),length(xval),j+length(xval)*(i-1))
                box on
               
                % Plot all data with negative A-values with filled markers
                if A_opt(s,r)<0
                    scatter(col(xval{j}(s,r)),col(Lomaxfunc{i}(LomaxParam{i}(s,r))),plottingsize(s,r),colors(s,:),markers{s,r},'filled');
                end 
                hold on
                
                % Plot all data with error bars
                if abs(LomaxParam{i}(s,r))<abs(LomaxErr{i}(s,r)) %if error crosses zero, make lower error bar black
                    errorbar(col(xval{j}(s,r)),col(Lomaxfunc{i}(LomaxParam{i}(s,r))),col(Lomaxfunc{i}(LomaxParam{i}(s,r)))-.001,...
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
        ylim(ylimit{i});
        eval(['set(gca,',axisscale{:},')']); 
        
        if j==1
            ylabel(LomaxParamname{i},'Interpreter','Latex','fontsize',14)
        end
        
        % Annotate with r-values calculated independently for A>0 (r+) and A<0 (r-)
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