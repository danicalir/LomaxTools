1. Download entire PNAS repository, which should contain:

	- createFit_LomaxR.m
	- distributions.m
	- make_hilofilteredrasters.m
	- make_rockdiststructure.m
	- make_rockstatstructure.m
	- plot_hillslopetransects.m
	- plot_rockdist_disentrainmentrate.m


2. Download the following scripts from Common repository and move into PNAS folder on your computer.

	- nangaussfilt.m
	- ndnanfilter.m
	- plothilofilter.m


3. Download all supplemental data files from SI Appendix. Unzip. Contents should include:
	
	Rasters:
	- H17c-avgz-1cm-bgplaneN.tif	
	- H20c-avgz-1cm-bgplaneN.tif	
     	- H28c-avgz-1cm-bgplaneN.tif	
	- N14c-avgz-1cm-bgplaneN.tif	
	- N20c-avgz-1cm-bgplaneN.tif	
	- N24c-avgz-1cm-bgplaneN.tif	
	- N39c-avgz-1cm-bgplaneN.tif	

	Experimental particle travel distances and right censor distances:
     	- HPB.rockdist.17deg.csv		
	- HPB.rockdist.20deg.csv		
	- HPB.rockdist.28deg.csv		
	- NOB.rockdist.14deg.csv		
	- NOB.rockdist.20deg.csv		
	- NOB.rockdist.24deg.csv		
	- NOB.rockdist.39deg.csv		
	- censorpoints.csv	*note: censor points have already been	

	Experimental particle characteristics:
	- rockstats.small.csv
	- rockstats.medium.csv
	- rockstats.large.csv


3. Create a directory within PNAS called ‘Data’ and move all supplemental data files into it.


4. Run the Matlab code from the PNAS working directory and in the following order:

	a. make_rockstats.m
		- makes "rockstats" structure containing experimental particle characteristics shown in SI Appendix, Table S2.

	b. make_filteredtopo.m
		- makes "topo" structure containing original rasters of surface-normal distance from bare-ground plane, filtered high-pass (short wavelength) and low-pass (long wavelength) rasters, and metadata.

	. make_rockdrops.m
		- makes "rockdrop" structure containing empirical particle travel distance distributions, raw experimental data and metadata
	

	. plot_rockdist_disentrainmentrate.m
		- plots

	. plot_hillslopetransects.m
		- plots all hillslope transects and particle travel distance histograms used in Fig. 2 and SI Appendix, Fig S2
		
	.	
		- 