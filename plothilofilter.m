function [h,a,s,c] = plothilofilter(raster,raster_hp,raster_lp,varargin)

% Generates figure showing unfiltered,  hi-pass and low-pass rasters.
% Optional subplot showing filter itself for scale.
% Outputs handle objects for plot surfaces, plotdim and contours.
%
% REQUIRED INPUTS:
%   raster = grid of original, unfiltered z values
%   raster_hp = grid of high-pass z values 
%   raster_lp = grid of low-pass z values corresponding with raster_hp
%
% OPTIONAL INPUTS:
%   x = 1D x-array
%   y = 1D y-array
%   dx = spacing in x
%   dy = spacing in y
%   filter = raster of filter used (will plot as separate subplot if entered) [default = []]
%   cmap = colormap array [64,3] (**default = 64-bit RedBlue)
%   contours = {N1,N2,N3}: number of contours on each raster figure
%           - N1 = contours on unfiltered surface figure
%           - N2 = contours on low-pass surface figure
%           - N3 = contours on high-pass surface figure
%           or {V1,V2,V3}: vectors for contour values on each figure
%           ** default (no value entered) = no contours [default = {}]
%   titletext = string with title of plot (shown on leftmost y-axis) [default = {}]
%   clims = limits for color map (use if all maps should have same colorbar) [default = {}]
%   plotdim = dimension (1=rows, 2=columns) along which to order subplots **default = columns (1x3 or 1x4)
% 
%
% OUTPUTS: 
%  Note: for all handle arrays, *(1)=unfiltered, *(2)=low-pass,  *(3)=high-pass, *(4)=filter
%   h = raster surface object handles
%   a = axis object handles
%   s = subplot object handles
%   c = contour object handles
%
% 
% EXAMPLE:
% [h,a,s,c] = plothilofilter(raster,raster_hp,raster_lp, 'filter',filter, 'x',x, 'y',y, 'contours',{5,-0.15:0.05:0.25,0}, 'titletext','Example Figure');
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% Written by Danica Roth, University of Oregon, October 2018.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% PARSE OPTIONAL INPUT ARGUMENTS
okargs={'filter' 'contours' 'titletext' 'clims' 'x' 'y' 'dx' 'dy' 'cmap' 'plotdim'};
defaults={[]   {}   {}   []   []   []   1   1   lbmap(64,'RedBlue')   2};
[filter,contours,titletext,clims,x,y,dx,dy,cmap,plotdim]=internal.stats.parseArgs(okargs, defaults, varargin{:});

minall = @(x) min(x,'all','omitnan');
if strcmpi(clims,'equal')
    clims=[min([raster,raster_hp,raster_lp],[],'all','omitnan') max([raster,raster_hp,raster_lp],[],'all','omitnan')];
end

if isempty(x)
    x=(0:size(raster,2)-1).*dx;
else
    dx=min(eldiff(x));
end
if isempty(y)
    y=(0:size(raster,1)-1).*dy;
else
    dy=min(eldiff(y));
end
[X,Y]=meshgrid(x,y);

%PLOT SURFACE
if plotdim==1
    s(1)=subplot(3+~isempty(filter),1,1);
elseif plotdim==2
    s(1)=subplot(1,3+~isempty(filter),1);
end
h(1)=surf(X,Y,raster);
a(1)=gca;
shading interp
view(0,90)
if ~isempty(clims)
    caxis(clims)
end
colormap(gca,cmap)
colorbar
hold on
if isempty(contours)==0
    [~,c(1)]=contour3(X,Y,raster,contours{1},'k','LineWidth',.05);
else
    c=[];
end
xlim([0 max(X,[],'all')]);
ylim([0 max(Y,[],'all')]);
title('Unfiltered')
if isempty(titletext)==0
    ylabel(titletext)
end

%PLOT LOW PASS SURFACE
if plotdim==1
    s(2)=subplot(3+~isempty(filter),1,2);
elseif plotdim==2
    s(2)=subplot(1,3+~isempty(filter),2);
end
h(2)=surf(X,Y,raster_lp);
a(2)=gca;
shading interp
view(0,90)
if ~isempty(clims)
    caxis(clims)
end
colormap(gca,cmap)
colorbar

hold on
if isempty(contours)==0
    [~,c(2)]=contour3(X,Y,raster_lp,contours{2},'k','LineWidth',.05);
end
xlim([0 max(X,[],'all')]);
ylim([0 max(Y,[],'all')]);
title('Low-Pass')

%PLOT HIGH PASS SURFACE
if plotdim==1
    s(3)=subplot(3+~isempty(filter),1,3);
elseif plotdim==2
    s(3)=subplot(1,3+~isempty(filter),3);
end
h(3)=surf(X,Y,raster_hp);
a(3)=gca;
shading interp
view(0,90)
if ~isempty(clims)
    caxis(clims)
end
colormap(gca,cmap)
colorbar();

hold on
if isempty(contours)==0
    [~,c(3)]=contour3(X,Y,raster_hp,contours{3},'k','LineWidth',.05);
end
xlim([0 max(X,[],'all')]);
ylim([0 max(Y,[],'all')]);
title('High-Pass')


%PLOT FILTER
if ~isempty(filter)
if plotdim==1
    s(4)=subplot(4,1,4);
elseif plotdim==2
    s(4)=subplot(1,4,4);
end
    xf=(0:size(filter,2)-1).*dx;
    yf=(0:size(filter,1)-1).*dy;
    [Xf,Yf]=meshgrid(xf,yf);
    h(4)=surf(Xf,Yf,filter);
    a(4)=gca;
    shading interp
    view(0,90)
    xlim([0 max(X,[],'all')]);
    ylim([0 max(Y,[],'all')]);
    colormap(gca,cmap)
    colorbar

    title('Filter')
end
