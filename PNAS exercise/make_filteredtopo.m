% Generate structure containing original rasters of surface-normal distance
% from bare-ground plane, filtered high-pass (short wavelength) and
% low-pass (long wavelength) rasters, roughness height distributions and
% metadata.  
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear

% Identify site names
topo.site={'B17' 'B20' 'B28' 'V14' 'V20' 'V24' 'V39'};

% Set gaussian filter parameters
filtsigma=100; %filter sigma value (in cm to match raster resolution)
hsize=1.*filtsigma; %filter window size

for i=1:length(topo.site)
    % Read in geotiff rasters and populate metadata fields in "topo" structure
    topo.geotiff_file=[topo.site{i},'c-avgz-1cm-bgplaneN.tif'];
    topo.slopedeg(i)=str2num(topo.site{i}(2:3));
    topo.slope(i)=tan(deg2rad(topo.slopedeg(i)));
    topo.raster{i}=imread(['Data/',topo.geotiff_file])';
    
    % Generate high- and low-pass rasters
    [topo.raster_lp{i},topo.raster_hp{i}] = nangaussfilt((topo.raster{i}),hsize,filtsigma,'plot','off'); 
    
    % Roughness calculations
    [topo.CDF_dhp{i},topo.dhp{i}]=ecdf(2.*abs(topo.raster_hp{i}(:))); %calculate empirical CDF of roughness heights d
    topo.d50(i,1:3)=prctile(2.*abs(topo.raster_hp{i}),50,'all'); %find median roughness height, d50
end

save('Data/topo.mat','topo','filtsigma','hsize');

disp('Topographic data saved to Data/topo.mat');