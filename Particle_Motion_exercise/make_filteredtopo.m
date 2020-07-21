% Generate structure containing original rasters of surface-normal distance
% from bare-ground plane, filtered high-pass (short wavelength) and
% low-pass (long wavelength) rasters, roughness height distributions and
% metadata.  
% 
% Note: If adjusting this code for different data, rasters (of surface-normal
% distances from bare-ground plane) for each site should be contained in
% geotiff files following naming format: L##.descriptors.tiff, where L is a
% site identifier (B or V) and ## is a 2-digit slope angle measured in
% degrees, and descriptors can be any additional naming details.
% All *.rockdist.csv files must have 3 columns, with each column containing
% data for the particle sizes described in the rockstats.mat data structure.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


clear

% Set gaussian filter parameters
filtsigma=100; %filter sigma value (in cm to match raster resolution)
hsize=1.*filtsigma; %filter window size

% Read in geotiff rasters and populate metadata fields in "topo" structure
files=dir('Data/*.tif');
topo.geotiff_file={files.name};
for i=1:length(topo.geotiff_file)
    eval(['topo.site{i}=''',topo.geotiff_file{i}(1:3),''';']);
    % Read in geotiff rasters and populate metadata fields in "topo" structure
    topo.slopedeg(i)=str2num(topo.site{i}(2:3));
    topo.slope(i)=tan(deg2rad(topo.slopedeg(i)));
    topo.raster{i}=imread(['Data/',topo.geotiff_file{i}]);
    [n,m]=size(topo.raster{i});
    if n>m
        topo.raster{i}=topo.raster{i}'; %reshape rasters so travel distance (assumed to be longer dimension) matches x dimension (columns)
    end
    
    % Generate high- and low-pass rasters (short/long wavelength topography) 
    [topo.raster_lp{i},topo.raster_hp{i}] = nangaussfilt((topo.raster{i}),hsize,filtsigma,'plot','off'); 
    
    % Roughness calculations
    [topo.CDF_dhp{i},topo.dhp{i}]=ecdf(2.*abs(topo.raster_hp{i}(:))); %calculate empirical CDF of roughness heights d
    topo.d50(i,1:3)=prctile(2.*abs(topo.raster_hp{i}),50,'all'); %find median roughness height, d50
end

save('Data/topo.mat','topo','filtsigma','hsize');

disp('Topographic data saved to Data/topo.mat');