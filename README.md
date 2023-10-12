# IPSL-PCM-Tools
Suite of tools to transform IPSL PCM outputs


###########################
########## VENUS ##########

createVenus_airglow.py myfile.nc
  --> create myfile_airglow3.nc containing airglow variables


createVenus_AS.py myfile.nc
  --> transform myfile.nc on a (lat ; lon) grid to myfile_AS.nc on a (angle to SS point ; azimuth to SS point) grid
  --> useful for suboslar-to-antisubsolar flow 


createVenus_ex.py myfile.nc
  --> create myfile_ex.nc with various useful variables (total wind speed, PV, pot temp, etc ...)

  
createVenus_supersonic.py  myfile.nc
  --> create myfile_ss.nc with various useful variables for supersonic flow (Mach number, speed of sound, eta ...)

  
createVenus_tt.py myfile.nc     
  --> transform myfile_tt.nc on time axis to myfile_tt.nc on a local time axis
  --> can give more information and is more vestaile than the IPSL Fortran tool localtime.F90

Some tools can be combined, e.g.
  > createVenus_supersonic.py myfile.nc
  > createVenus_AS.py myfile_ss.nc
will create myfile_ss_AS.nc

