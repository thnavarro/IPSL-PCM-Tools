#!/usr/bin/env python


from    optparse   import  OptionParser
import  numpy      as      np
from    netCDF4    import  Dataset
import  gc

# T Navarro - 2022

# create_tt: create a separate file of Venus (LMD) GCM outputs on a LT grid
# (time; lat ; lon) --> (LT; lat; lon)


parser = OptionParser()

(opt,args) = parser.parse_args()

if args is None:
    print('Stop here! I need a file as argument(s)!'); exit()
else:
    fname = args[0]

nc = Dataset(fname,'r+')
ncv = nc.variables


if '_tt.nc' in fname:
  print('I think this is already a tt file')
  print('Abort!!')
  exit()

if '_A' in fname:
  namelon = 'longitude'
  namelat = 'latitude'
  nametime = 'Time'
#  namealt = 'altitude'
else:
  namelon = 'lon'
  namelat = 'lat'
  nametime = 'time_counter'
#  namealt = 'presnivs'


nc = Dataset(fname,'r+')
ncv = nc.variables


lat  = ncv[namelat][:]
lon  = ncv[namelon][:]
#alt  = ncv[namealt][:]
time  = ncv[nametime][:]

time2 = (time-2*time[0]+time[1])/10087000.


#aslon = -(time2%1)*360.
#aslon[aslon<180] = aslon[aslon<180] + 360.

#lon2 = lon*24./360 + 12.
lon2 = lon + 0.
lat2 = lat + 0.


#lt = np.linspace(0.,24.*int(time2[-1]),num=len(time))
lt = np.linspace(0.,24.,num=len(time2[time2<=1]))

# for percentile in amplitude caclulation
bins = [80,85,90,95,99]

fnameout=fname[0:len(fname)-3]+'_tt.nc' # tt is for thermal tide

nc2 = Dataset(fnameout, 'w', format='NETCDF3_64BIT')
nc2.createDimension(namelon,len(lon2))
nc2.createDimension(namelat,len(lat2))
#nc2.createDimension(namealt,len(alt))
nc2.createDimension(nametime,len(time2))
nc2.createDimension('loctime',len(lt))
nc2.createDimension('bins',len(bins))

lon_ = nc2.createVariable(namelon, 'f', (namelon,))
lon_[:] = lon2[:]
lat_ = nc2.createVariable(namelat, 'f', (namelat,))
lat_[:] = lat2[:]
#alt_ = nc2.createVariable(namealt, 'f', (namealt,))
#alt_[:] = alt[:]
time_ = nc2.createVariable(nametime, 'f', (nametime,))
time_[:] = time2[:]
lt_ = nc2.createVariable('loctime', 'f', ('loctime',))
lt_[:] = lt[:]
bins_ = nc2.createVariable('bins', 'f', ('bins',))
bins_[:] = bins[:]


zevars = [
#'tops',
'psol',
#'temp',
'speed',
#'realw',
#'co2'
#'o2',
#'o',
#'o2dg'
#'bv2',
#'pressure',
#'height',
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
          print('not yet')
          exit()
#         field = ncv[zevar][:,:,:,:]
#         field2 = np.zeros((len(time2),len(alt),len(lat2),len(lon2)))
     elif ndim == 3:
         field     = ncv[zevar][:,:,:]
         field_lt  = np.zeros((len(time2),len(lat2),len(lon2)))
         field_ltm = np.zeros((len(lt),len(lat2),len(lon2)))
         field_lts = np.zeros((len(lt),len(lat2),len(lon2)))
         field_min = np.zeros((len(lt),len(lat2),len(lon2)))
         field_max = np.zeros((len(lt),len(lat2),len(lon2)))
         for i in range(len(lon2)):
            # temps pour lequel LT = minuit 
            midnight = (time2 - lon[i]/360.)%1
            midnight0 = np.argmin(abs(midnight))
           # print(midnight0)
            tmp  = np.roll(field[:,:,i],-midnight0,axis=0)
            field_lt[:,:,i] = tmp
         sigma = time2[1] -time2[0]
         anom = field_lt - field_lt.mean(0)
         for i in range(len(time2[time2<=1])):
            inds = abs(time2%1 - time2[i]%1)<sigma/3.
            field_ltm[i,:,:] = anom[inds,:,:].mean(0)
            field_lts[i,:,:] = anom[inds,:,:].std(0)
            field_min[i,:,:] = anom[inds,:,:].min(0)
            field_max[i,:,:] = anom[inds,:,:].max(0)
         amp = np.zeros((len(bins),len(lat2),len(lon2)))
         for i in range(len(bins)):
            amp[i,:,:] = np.percentile(field,bins[i],axis=0) - np.percentile(field,100-bins[i],axis=0)
         amp_max = field_max.max(0) - field_min.min(0)
         amp_min = field_min.max(0) - field_max.min(0)
         amp_ave = field_ltm.max(0) - field_ltm.min(0)
         amp_amp = amp_max - amp_min 
         amp_3sig = 3*field.std(0) 

 
     if   ndim == 4:
       #field_ = nc2.createVariable(zevar,'f', (nametime,namealt,namelat,namelon,))
       print('no code here yet')
     elif ndim == 3:
       field_  = nc2.createVariable(zevar,'f', (nametime,namelat,namelon,))
       field_m = nc2.createVariable(zevar+'_m','f', ('loctime',namelat,namelon,))
       field_s = nc2.createVariable(zevar+'_s','f', ('loctime',namelat,namelon,))
       field_a = nc2.createVariable(zevar+'_a','f', ('bins',namelat,namelon,))
       field_amax = nc2.createVariable(zevar+'_amax','f', (namelat,namelon,))
       field_amin = nc2.createVariable(zevar+'_amin','f', (namelat,namelon,))
       field_aave = nc2.createVariable(zevar+'_aave','f', (namelat,namelon,))
       field_aamp = nc2.createVariable(zevar+'_aamp','f', (namelat,namelon,))
       field_3sig = nc2.createVariable(zevar+'_3sig','f', (namelat,namelon,))

     field_[:] = field_lt[:]
     field_.info = 'Field created by script createVenus_tt.py'
     field_m[:] = field_ltm[:]
     field_m.info = 'Field created by script createVenus_tt.py'
     field_s[:] = field_lts[:]
     field_s.info = 'Field created by script createVenus_tt.py'

     field_a[:] = amp[:]
     field_a.info = 'Field created by script createVenus_tt.py'

     field_amax[:] = amp_max[:]
     field_amax.info = 'Field created by script createVenus_tt.py'

     field_amin[:] = amp_min[:]
     field_amin.info = 'Field created by script createVenus_tt.py'

     field_aave[:] = amp_ave[:]
     field_aave.info = 'Field created by script createVenus_tt.py'

     field_aamp[:] = amp_amp[:]
     field_aamp.info = 'Field created by script createVenus_tt.py'

     field_3sig[:] = amp_3sig[:]
     field_3sig.info = 'Field created by script createVenus_tt.py'

     del field_lt,field,field_
     gc.collect()


nc.close()
nc2.close()



