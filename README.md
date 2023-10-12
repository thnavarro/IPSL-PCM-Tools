# IPSL-PCM-Tools
Suite of tools to transform IPSL PCM outputs


## VENUS  

Tip: some tools can be combined, e.g.,  
  > createVenus_supersonic.py myfile.nc  
  > createVenus_AS.py myfile_ss.nc

will create myfile_ss_AS.nc  

### createVenus_airglow.py myfile.nc
Create myfile_airglow3.nc containing airglow variables.

### createVenus_AS.py myfile.nc  
Transform myfile.nc on a (lat ; lon) grid to myfile_AS.nc on a (angle to SS point ; azimuth to SS point) grid.  
Useful for suboslar-to-antisubsolar flow. 

### createVenus_ex.py myfile.nc  
Create myfile_ex.nc with various useful variables (total wind speed, PV, pot temp, etc ...)  

### createVenus_supersonic.py  myfile.nc  
Create myfile_ss.nc with various useful variables for supersonic flow (Mach number, speed of sound, eta ...)  
   
### CreateVenus_tt.py myfile.nc  
Transform myfile.nc on time axis to myfile_tt.nc on a local time axis.  
It can give more information and is more vestaile than the IPSL Fortran tool localtime.F90  
  
  


