# LomaxTools

These Matlab functions perform statistical and spectral analysis on empirical observations and gridded topographic data as described in Roth et al. (submitted to PNAS, 2020): Particle Motion on Burned and Vegetated Hillslopes<sup><b><i>(1)</i></b></sup>. 

The scripts contained in the Particle_Motion_exercise directory demonstrate the implementation of these functions, as well as initial data file preparation and figure generation. Data from this publication are contained in https://github.com/danicalir/Data_ParticleMotion. Please acknowledge the use of this software in any publications by citing this paper and software release.

[![DOI](https://zenodo.org/badge/281159851.svg)](https://zenodo.org/badge/latestdoi/281159851)



To run the example in ./Particle_Motion_exercise:
1) Download data files from https://github.com/danicalir/Data_ParticleMotion and move or copy the entire 'Data' directory into the 'Particle_Motion_exercise' directory.
2) Open and run run_allcode.mlx in Matlab from the 'Particle_Motion_exercise' working directory.

See comments in individual Matlab scripts for more details. For help with any function, type 'help [function name]' at the Matlab prompt.

All code copyright (C) 2020 Danica Roth (droth@mines.edu), except where otherwise indicated. These programs are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation. You should have received a copy of the GNU General Public License along with these programs. If not, see http://www.gnu.org/licenses.


<b>Included files:</b>

<i>Matlab functions:</i>

- createfit.m

- distributions.m

- fitgpdistributions.m

- lomaxbootstrap.m 

- lomaxopt.m

- nangaussfilt.m

- ndnanfilter.m<sup><b><i>(2)</i></b></sup>

- plothilofilter.m

- ssepareto.m


<br>


<i>Matlab scripts contained in ./Particle_Motion_exercise:</i>

- make_filteredtopo.m

- make_rockdrops.m

- make_rockstats.m

- plot_allfigures.m

- plot_comparegpfits.m

- run_allcode.mlx
  
- run_lomaxbootstraperr.m
  
- run_lomaxoptimization.m


<br>
<br>

<b>Notes:</b>

<b>(1)</b> This project was supported by funding from the National Science Foundation (EAR-1420898 and Postdoctoral Fellowship 1625311).

<b>(2)</b> Written by Carlos Adrian Vargas Aguilera. Obtained from the Matlab Central File Exchange:
<br>
&nbsp;&nbsp;&nbsp;&nbsp;https://www.mathworks.com/matlabcentral/fileexchange/20417-ndnanfilter-m
<br>
&nbsp;&nbsp;&nbsp;&nbsp;Distributed under the BSD license (see m-file header for more information).
