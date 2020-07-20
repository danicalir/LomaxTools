clear
close all
load('Data/results.mat')
load('Data/topo.mat');

name = {'original' 'low-pass' 'high-pass'};
iorder=[2 4 6 1 3 5 7]; %reorder X, R, Censorpoint from figure order from plot_rockdst_disentrainmentrate_FINAL.m
colors=[237 177 32; 179 177 32; 119 172 48; 119 172 48; 77 190 238; 0 131 205; 126 47 142]./255;%colors for plotting hi-pass
colors=colors(iorder,:);%reorder (originally from fig.2 order)
markersize=[1.7 4.5 7.2].^1.3;
fc=[     0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];%colors for histograms  
yt=[300;...%transect y(element) number (manually selected for each site)
    150;...
    120;...
    94;...
    220;...
    155;...
    304];
isite='BBBNNNN';
hnum=[10 30 30]; %number of histogram bins for different particle sizes
theta=0.1; %theta (left truncation) at all sites = 0.1 m

CensorPoint=cell2mat(cellfun(@(v)v.distributions.CensorPoint,results,'Uniform',0));
for i = 1:7
    nx=size(topo.raster{i},2); 
    ny=size(topo.raster{i},1);
    y = [1:ny]./100;   
    x = [1:nx]./100;
    z = topo.raster{i}(yt(i),:);
    zs = z-(topo.slope(i).*x); %re-trend elevations in horizontal reference frame
    x = x./cosd(topo.slopedeg(i))-theta; %convert x to along-ground reference frame (matching travel distances) and shift by truncation distance (theta)
    z0 = interp1(x(~isnan(zs)),zs(~isnan(zs)),0,'linear','extrap'); %extrapolate linearly to find z(x=0) to rezero z (for sites with z(x=0)=NaN)
    zs = zs-z0; %rezero z to x=0 (theta-shifted starting line) elevation
    cp(i)=max(CensorPoint(i,:)-theta); %make i-specific cp and shift by theta for easy reference later
   
    %FILTER RASTER
    raster={topo.raster{i} topo.raster_lp{i} topo.raster_hp{i}};
    
    %Plot figures
    figure(i+10)
    clf
    
    %Plot the topographic transects 
    s1 = subplot(2,1,1);
    yyaxis left
    plot(x,zs,'k','linewidth',[2]); hold on
    ylabel('elevation below starting point [m]')
    
    %Set lower z-limit for plotting
    if i==1
        zlimval = [-11 0]; %manually set first site because raster does not extend to max travel distance
    else
        zlimval = [interp1(x(~isnan(zs)),zs(~isnan(zs)),1.05*min(max(x),cp(i)),'linear','extrap') 0];
    end
    
    %Plot the hi-pass topography
    yyaxis right
    plot(x,raster{2}(yt(i),:),'-s','linewidth',5,'Color',colors(i,:));
    plot(x,raster{3}(yt(i),:),'linewidth',[1],'Color',colors(i,:));
    ylabel('filtered topography [m]')    
    set(s1.YAxis(1),'Color','k','Limits',zlimval)
    set(s1.YAxis(2),'Color',colors(i,:),'Limits',[-.75 .75])
    set(s1.XAxis,'Limits',[0 1.05*cp(i)]);
    title(topo.site{i}(1:3));
    
    %Plot close-up inset
    axes('Position',[.64 .84 .25 .07])
    box on
    plot(x,raster{3}(yt(i),:),'linewidth',[1],'Color',colors(i,:));
    ylabel('d [m]')
    xlabel('x [m]')
    ylim([-.3 .3])
    xlim([0 2.3])
    ax2=gca;
    ax2.YAxis(1).Color = colors(i,:);

    %Plot histograms of travel distances
    subplot(2,1,2)
    for r=[3 2 1]
        xdata{r} = results{i,r}.data.x(results{i,r}.data.x>0.1)-0.1;
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


% % % Position final figures for viewing
% a = [-231;
%      -1;
%      229;
%      461;
%      692;
%      923;
%     1156];
% b = 1051;
% c = 638;
% d = 984;
% 
% for i=fliplr(1:7)
%     figure(i+10)
%     set(gcf,'position',[a(i) b c d]);
% end
