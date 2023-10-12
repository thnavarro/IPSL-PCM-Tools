#!/usr/bin/env python


from    optparse   import  OptionParser
import  numpy      as      np
from    netCDF4    import  Dataset
import  gc

# T Navarro - 2020

# create_AS: create a separate file of Venus (LMD) GCM outputs on a grid to understand the AntiSolar (AS) point
# (lat ; lon) --> (angle to SS point ; azimuth to SS point)
# see Koll et al. for grid transformation of tidally locked exoplanets


parser = OptionParser()

(opt,args) = parser.parse_args()

if args is None:
    print('Stop here! I need a file as argument(s)!'); exit()
else:
    fname = args[0]

nc = Dataset(fname,'r+')
ncv = nc.variables




if '_A' in fname:
  namelon = 'longitude'
  namelat = 'latitude'
  nametime = 'Time'
  namealt = 'altitude'
else:
  namelon = 'lon'
  namelat = 'lat'
  nametime = 'time_counter'
  namealt = 'presnivs'


nc = Dataset(fname,'r+')
ncv = nc.variables


lat  = ncv[namelat][:]
lon  = ncv[namelon][:]
alt  = ncv[namealt][:]
time  = ncv[nametime][:]

time2 = (time-2*time[0]+time[1])/10087000.


aslon = -time2*360.
#aslon[aslon<180] = aslon[aslon<180] + 360.


#lon2 = lon + 0.
#lat2 = lat + 0.
lon2 = lon[::2] + 0.
lat2 = lat[::2] + 0.


fnameout=fname[0:len(fname)-3]+'_AS.nc' # AS is for AS point

nc2 = Dataset(fnameout, 'w', format='NETCDF3_64BIT')
nc2.createDimension(namelon,len(lon2))
nc2.createDimension(namelat,len(lat2))
nc2.createDimension(namealt,len(alt))
nc2.createDimension(nametime,len(time))

lon_ = nc2.createVariable(namelon, 'f', (namelon,))
lon_[:] = lon2[:]
lat_ = nc2.createVariable(namelat, 'f', (namelat,))
lat_[:] = lat2[:]
alt_ = nc2.createVariable(namealt, 'f', (namealt,))
alt_[:] = alt[:]
time_ = nc2.createVariable(nametime, 'f', (nametime,))
time_[:] = time[:]


zevars = [
'tops',
'psol',
'temp',
'speed',
#'realw',
#'co2'
#'o2',
#'o',
#'o2dg'
#'bv2',
'pressure',
'height',
'ps',
'mach',
'cs',
'eta',
]
# note : vitu and vitv are grid-dependent, would require a transformation.


for zevar in zevars:
  if zevar in ncv:

     ndim = len(ncv[zevar].shape)
     if   ndim == 4:
         field = ncv[zevar][:,:,:,:]
         field2 = np.zeros((len(time2),len(alt),len(lat2),len(lon2)))
     elif ndim == 3:
         field = ncv[zevar][:,:,:]
         field2 = np.zeros((len(time2),len(lat2),len(lon2)))
     #field2 = np.zeros_like(field)
   
     for t in range(len(time2)):
       for i in range(len(lat2)):
         for j in range(len(lon2)):
            lat_t = np.arcsin(np.cos(lat2[i]*np.pi/180.)*np.cos(lon2[j]*np.pi/180.))
            lon_t = np.arctan2(np.sin(lon2[j]*np.pi/180.),np.tan(lat2[i]*np.pi/180.)) - aslon[t]*np.pi/180.
  
            lat_t = lat_t%(2*np.pi) 
            if lat_t > np.pi : lat_t = lat_t-2*np.pi 
            lon_t = lon_t%(2*np.pi) 
            if lon_t > np.pi : lon_t = lon_t-2*np.pi 
           
            ilat = np.argmin(abs(lat*np.pi/180.-lat_t))
            ilon = np.argmin(abs(lon*np.pi/180.-lon_t))
   
           # if i == 0 and j == 0:
           #   print(t,i,j)
           #   print(ilat,ilon)
           #   print(lat2[i],lon2[j])
           #   print(lat_t*180/np.pi,lon_t*180/np.pi,aslon[t])
           #   print('')
   
            if   ndim == 4: field2[t,:,i,j] = field[t,:,ilat,ilon]
            elif ndim == 3: field2[t,i,j]   = field[t,ilat,ilon] 
   
     if   ndim == 4:
       field_ = nc2.createVariable(zevar,'f', (nametime,namealt,namelat,namelon,))
     elif ndim == 3:
       field_ = nc2.createVariable(zevar,'f', (nametime,namelat,namelon,))
     field_[:] = field2[:]
     field_.info = 'Field created by script createVenus_AS.py'

     del field2,field,field_
     gc.collect()


nc.close()
nc2.close()



