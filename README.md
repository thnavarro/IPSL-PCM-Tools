# IPSL-PCM-Tools
Suite of tools to transform IPSL PCM outputs


## VENUS  

Syntax is:  
> script.py myfile.nc

with myfile.nc a IPSL PCM histmth output, that will create a new file myfile_XX.nc  
  
Tip: some tools can be combined, e.g.,  
  > createVenus_supersonic.py myfile.nc  
  > createVenus_AS.py myfile_ss.nc

will create myfile_ss_AS.nc  

### createVenus_airglow.py
Create myfile_airglow3.nc containing airglow variables.

### createVenus_AS.py
Transform myfile.nc on a (lat ; lon) grid to myfile_AS.nc on a (angle to SS point ; azimuth to SS point) grid.  
Useful for subsolar-to-antisubsolar flow. 

### createVenus_ex.py
Create myfile_ex.nc with various useful variables (total wind speed, PV, pot temp, etc ...)  

### createVenus_supersonic.py
Create myfile_ss.nc with various useful variables for supersonic flow (Mach number, speed of sound, eta ...)  
   
### createVenus_tt.py  
Transform myfile_tt.nc on time axis to myfile_tt.nc on a local time axis.  
It can give more information and is more vestaile than the IPSL Fortran tool localtime.F90  
  
  


