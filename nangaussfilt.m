function [raster_lp,raster_hp,h,a,s,c] = nangaussfilt(raster,hsize,sigma,varargin)

%Filters input raster with gaussian, optionally plots filter, unfiltered raster,
%and resulting hi-pass and low-pass rasters
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   %
% REQUIRED INPUTS:
%   raster = raster grid of z values 
%   hsize = filter size (in elements)
%   sigma = gaussian sigma value (in elements)
%
% OPTIONAL INPUTS:
%   plot = 'on' to generate figure [default = 'off']
%
%  Note: Remaining optional inputs are used for plothilofilter function and only
%  needed if plotting
%       x = 1D x-array
%       y = 1D y-array
%       dx = spacing in x
%       dy = spacing in y
%       cmap = colormap array [64,3] (**default = 64-bit RedBlue)
%       contours = {N1,N2,N3}: number of contours on each raster figure
%               - N1 = contours on unfiltered surface figure
%               - N2 = contours on low-pass surface figure
%               - N3 = contours on high-pass surface figure
%               or {V1,V2,V3}: vectors for contour values on each figure
%               ** default (no value entered) = no contours [default = {}]
%       titletext = string with title of plot (shown on leftmost y-axis) [default = {}]
%       clims = limits for color map (use if all maps should have same colorbar) [default = {}]
%       plotdim = dimension (1=rows, 2=columns) along which to order subplots **default = columns (1x3 or 1x4)
% 
%     
% OUTPUTS:
%   raster_lp = low pass raster
%   raster_h = high-pass raster
%  Note: for all handle arrays, *(1)=unfiltered (left), *(2)=lo-pass (middle), *(3)=hi-pass (right)
%   h = raster surface object handles
%   s = subplot handles
%   c = contour object handles
% 
% EXAMPLE:
% raster = rasterobjectname;
% hsize = z,hsize,filtsigma,'plot','on');
% [raster_lp,raster_hp,h,s,c] = nangaussfilt(zraster,10,10, 'filter',filter, 'contours',{5,-0.15:0.05:0.25,0}, 'titletext','Example Figure');
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   
% Written by Danica Roth, University of Oregon, October 2018.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
 
  
%% PARSE OPTIONAL INPUT ARGUMENTS
okargs={'plot' 'contours' 'titletext' 'clims' 'x' 'y' 'dx' 'dy' 'cmap' 'plotdim'};
defaults={'off'   {}   {}   []   []   []   []   []   []   []};
[plotstat,contours,titletext,clims,x,y,dx,dy,cmap]=internal.stats.parseArgs(okargs, defaults, varargin{:});

%% FILTER
if mod(hsize,2)==0
    hsize=hsize+1; %make odd filter size to ensure existence of center element
end
filter=fspecial('gaussian',hsize,sigma); %create gaussian filter
F=(hsize-1)/2; %create sliding window size over which to apply filter
[raster_lp]=ndnanfilter(double(raster),filter,[F,F]); %low pass surface
raster_hp=double(raster)-double(raster_lp); %high pass surface

%% PLOT FILTER (if 'plot' set to 'on')
if strcmpi(plotstat,'on')
    varargin(find(strcmp(varargin,'plot' )):find(strcmp(varargin,'plot' ))+1)=[]; %delete plot arguments from varargin so it can be passed to plothilofilter
    [h,a,s,c] = plothilofilter(raster,raster_hp,raster_lp, 'filter',filter, varargin{:});%'contours',contours, 'titletext',titletext,'clims',clims); 
else
    h=[];
    a=[];
    s=[];
    c=[];
end
