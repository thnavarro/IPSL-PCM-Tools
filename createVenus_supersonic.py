#!/usr/bin/env python


from    optparse   import  OptionParser
import  numpy      as      np
from    netCDF4    import  Dataset
import  gc

# T Navarro - 2020

# createVenus_supersonic:
# create an extended file of a Venus (LMD) GCM output.
# the separate file contains new variables relevant to supersonic flow (mach, etc ...)


######################
############ Constants
nu = 0.35
cp0 = 1000.
pref = 9.2e6
temp0 = 460. 
R = 8.314
g = 8.87
a = 6050e3
gamma = 1.3
######################
kappa0 = R/cp0 
######################

#################################################################################
#################################################################################
ap = np.array([0, 656.152424151957, 2610.736779638, 7673.04963622473,  \
    17327.74808857, 32955.9002381228, 55690.3398135276, 86255.7068921984, \
    124733.193259007, 170330.325009327, 221231.5654731, 274615.023359711, \
    326923.721714309, 374262.972598652, 412921.789530273, 439755.923293322, \
    451966.777242386, 450005.898678512, 434661.362963411, 401194.186762591, \
    328917.667350411, 237439.767615529, 170133.103413588, 121905.717423763, \
    87349.2912061494, 62588.5220813229, 44846.6499929582, 32134.0391501656, \
    23025.0524969264, 16498.1763071941, 11821.4636448756, 8470.45154178382, \
    6069.3456714494, 4348.87759627234, 3116.10732860692, 2232.78855630705, \
    1599.86277733637, 1146.35150211961, 821.396680616519, 588.556298745786, \
    421.718975708101, 302.174775044844, 216.51766937686, 155.141650733953, \
    111.163845616808, 79.6523574513065, 51.4968533417143, 28.2620718691288, \
    14.2948315243741, 6.1098240738474, \
    2.94853179385991, 1.79778664326264, \
    1.09041271876065, 0.661368745668942, 0.4011404216239, 0.243303964564948, \
    0.147571314138277, 0.0895065265189495, 0.0542884525781248, \
    0.032927610956988, 0.0199716055965029, 0.0121133911179674, \
    0.00734714310613793, 0.00445626755516896, 0.00270286290009264, \
    0.00163936921790599, 0.000994327693249103, 0.000603090231756919, \
    0.000365792716133769, 0.000221864497434691, 0.000134567619995875, \
    8.16193873320571e-05, 4.95046608438535e-05, 3.00260946004727e-05, \
    1.82117469666186e-05, 1.10459829021827e-05, 6.69972729683537e-06, \
    4.06359001724429e-06, 0.])

bp = np.array([1, 0.998688124935515, 0.994786750274124, 0.984711592615347, \
    0.965631408487589, 0.935085441244573, 0.891310686248881, \
    0.833485848239643, 0.76198545182255, 0.678516347149576, \
    0.586092910995097, 0.48881848341266, 0.391400247225106, \
    0.298631253175517, 0.214741704129835, 0.142900858694672, \
    0.0875163783425901, 0.048995779821085, 0.0229092855468057, \
    0.00666028808697784, 0.000266994865634196, 5.38254462282418e-08, \
    2.69344453442487e-15, 1.62378981504663e-29, 3.26745318386404e-57, \
    3.70832564210981e-111, 3.15260821245907e-216, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

#################################################################################


parser = OptionParser()

(opt,args) = parser.parse_args()

if len(args)==0:
    print('Stop here! I need a file as argument(s)!'); exit()
else:
    fname = args[0]

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

#it = np.arange(5)
it = np.arange(len(time))
time = time[it]

vitu = ncv['vitu'][it,:,:,:]
vitv = ncv['vitv'][it,:,:,:]
temp = ncv['temp'][it,:,:,:]
mmean = ncv['mmean'][it,:,:,:]*1e-3

nc.close()

speed = np.sqrt(vitu*vitu + vitv*vitv)


mlon,mlat = np.meshgrid(lon,lat)
mlat = np.cos(mlat*np.pi/180.)


dlon = (lon[1]-lon[0])
dlat = (lat[1]-lat[0])


#########################################
######### Compute divergence

eta = np.zeros_like(vitu)
#dv = (vitv[:,:,2:,1:-1]*mlat[np.newaxis,np.newaxis,2:,1:-1] - vitv[:,:,:-2,1:-1]*mlat[np.newaxis,np.newaxis,2:,1:-1])
dv = (vitv[:,:,2:,1:-1] - vitv[:,:,:-2,1:-1])
du = (vitu[:,:,1:-1,2:] - vitu[:,:,1:-1,:-2])
eta[:,:,1:-1,1:-1] = (du+dv) #/(a*mlat[np.newaxis,np.newaxis,1:-1,1:-1])
eta[:,:,1,:] = 0.
eta[:,:,-2,:] = 0.

#debug
#du2 = np.zeros_like(vitu)
#dv2 = np.zeros_like(vitu)
#du2[:,:,1:-1,1:-1] = du
#dv2[:,:,1:-1,1:-1] = dv


cs = np.sqrt(gamma*R*temp/mmean)
eta = -eta/cs 

mach = speed/cs


del temp,vitu,vitv,mmean #,cs
gc.collect()




fnameout=fname[0:len(fname)-3]+'_ss.nc' # ss is for supersonic

nc2 = Dataset(fnameout, 'w', format='NETCDF3_64BIT')
nc2.createDimension(namelon,len(lon))
nc2.createDimension(namelat,len(lat))
nc2.createDimension(namealt,len(alt))
nc2.createDimension(nametime,len(time))

lon_ = nc2.createVariable(namelon, 'f', (namelon,))
lon_[:] = lon[:]
lat_ = nc2.createVariable(namelat, 'f', (namelat,))
lat_[:] = lat[:]
alt_ = nc2.createVariable(namealt, 'f', (namealt,))
alt_[:] = alt[:]
time_ = nc2.createVariable(nametime, 'f', (nametime,))
time_[:] = time[:]


eta_   = nc2.createVariable('eta','f',  (nametime,namealt,namelat,namelon,))
#du_   = nc2.createVariable('du','f',  (nametime,namealt,namelat,namelon,))
#dv_   = nc2.createVariable('dv','f',  (nametime,namealt,namelat,namelon,))
cs_   = nc2.createVariable('cs','f',  (nametime,namealt,namelat,namelon,))
mach_   = nc2.createVariable('mach','f',  (nametime,namealt,namelat,namelon,))

###

eta_[:] = eta[:]
eta_.info = 'eta diagnostic'
#du_[:] = du2[:]
#dv_[:] = dv2[:]
cs_[:] = cs[:]
mach_[:] = mach[:]

nc2.close()



