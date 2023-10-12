#!/usr/bin/env python


from    optparse   import  OptionParser
import  numpy      as      np
from    netCDF4    import  Dataset
import  gc

# T Navarro - 2020

# createVenus_ex:
# create an extended file of a Venus (LMD) GCM output.
# the separate file contains new variables (wind speed, PV, pot tempe, etc ...)


######################
############ Constants
nu = 0.35
cp0 = 1000.
pref = 9.2e6
temp0 = 460. 
R = 191.4 # oui mais non car la composition change
g = 8.87
a = 6050e3
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

vitu = ncv['vitu'][:,:,:,:]
vitv = ncv['vitv'][:,:,:,:]
temp = ncv['temp'][:,:,:,:]
vitw = ncv['vitw'][:,:,:,:]


speed = np.sqrt(vitu*vitu + vitv*vitv)


mlon,mlat = np.meshgrid(lon,lat)
mlat = np.cos(mlat*np.pi/180.)


dlon = (lon[1]-lon[0])
dlat = (lat[1]-lat[0])

#########################################
######### Compute vorticity

rot = np.zeros_like(vitu)
du = (vitu[:,:,2:,1:-1]*mlat[np.newaxis,np.newaxis,2:,1:-1] - vitu[:,:,:-2,1:-1]*mlat[np.newaxis,np.newaxis,2:,1:-1])/dlat
dv = (vitv[:,:,1:-1,2:] - vitv[:,:,1:-1,:-2])/dlon
rot[:,:,1:-1,1:-1] = (-du-dv)/(a*mlat[np.newaxis,np.newaxis,1:-1,1:-1])
rot[:,:,1,:] = 0.
rot[:,:,-2,:] = 0.


#########################################
######### Compute divergence

# div = np.zeros_like(vitu)
# dv = (vitv[:,:,2:,1:-1]*mlat[np.newaxis,np.newaxis,2:,1:-1] - vitv[:,:,:-2,1:-1]*mlat[np.newaxis,np.newaxis,2:,1:-1])/dlat
# du = (vitu[:,:,1:-1,2:] - vitu[:,:,1:-1,:-2])/dlon
# div[:,:,1:-1,1:-1] = (du-dv)/(a*mlat[np.newaxis,np.newaxis,1:-1,1:-1])
# div[:,:,1,:] = 0.
# div[:,:,-2,:] = 0.



#########################################
######### Compute teta

if 'pres' in ncv:
  pres = ncv['pres'][:,:,:,:]
  ps = ncv['psol'][:,:,:]
  pres_interlay = np.zeros((len(time),len(ap),len(lat),len(lon)))
  pres_interlay[:,1:-1,:,:] = 0.5*(pres[:,1:,:,:] + pres[:,:-1,:,:])
  pres_interlay[:,0,:,:] = ps
  alts = -R*temp*np.diff(np.log(pres_interlay[:,:,:,:]),axis=1)/g
  alts = np.cumsum(alts,axis=1)
elif 'pressure' in ncv:
  pres = ncv['pressure'][:,:,:,:]
else:
  ps = ncv['psol'][:,:,:]

  pstiled = ps[:,np.newaxis,:,:]
  aptiled = ap[np.newaxis,:,np.newaxis,np.newaxis]
  bptiled = bp[np.newaxis,:,np.newaxis,np.newaxis]

  pres_interlay = pstiled*bptiled + aptiled
  gc.collect()
  pres = 0.5*(pres_interlay[:,1:len(ap),:,:] + pres_interlay[:,0:len(ap)-1,:,:])
  alts = -R*temp*np.diff(np.log(pres_interlay[:,:,:,:]),axis=1)/g
  alts = np.cumsum(alts,axis=1)
  del pres_interlay,pstiled,aptiled,bptiled
  gc.collect()
  print('createVenus_ex.py: 3D pressure reconstructed from surface pressure.')
  #exit()

teta = np.power(temp,nu) + nu*np.power(temp0,nu)*kappa0*np.log(pref/pres)
teta = np.power(teta,1./nu)

#dteta = np.zeros_like(teta)
#dteta[:,1:,:,:] = (teta[:,:-1,:,:] - teta[:,1:,:,:])/(pres[:,1:,:,:] - pres[:,:-1,:,:])


#########################################
######### Compute Ertel PV

#pv = rot*dteta


realw = vitw*R*temp/(g*pres)

nc.close()


###########################################
######### Compute Brunt Vaisala Frequency


#bv2 = np.zeros_like(vitu)
#bv2[:,1:nlev,:,:] = 2*g*(-teta[:,1:nlev,:,:]+teta[:,0:nlev-1,:,:])  \
#                      /  (teta[:,1:nlev,:,:]+teta[:,0:nlev-1,:,:]) \
#                      /  (pres[:,1:nlev,:,:]-pres[:,0:nlev-1,:,:])
#bv2 = g*bv2*dpdz

#bv[bv2>=0] = np.sqrt(bv[bv2>=0])
#bv[bv2<0] = -np.sqrt(-bv[bv2<0])


#del pres,dteta,temp,vitu,vitv,vitw
del pres,temp,vitu,vitv,vitw #,bv2
gc.collect()




fnameout=fname[0:len(fname)-3]+'_ex.nc' # ex is for extended

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


alts_ = nc2.createVariable('height','f', (nametime,namealt,namelat,namelon,))
speed_ = nc2.createVariable('speed','f', (nametime,namealt,namelat,namelon,))
rot_   = nc2.createVariable('vort','f',  (nametime,namealt,namelat,namelon,))
div_   = nc2.createVariable('div','f',  (nametime,namealt,namelat,namelon,))
teta_  = nc2.createVariable('teta','f',  (nametime,namealt,namelat,namelon,))
#pv_    = nc2.createVariable('pv','f',    (nametime,namealt,namelat,namelon,))
realw_ = nc2.createVariable('realw','f', (nametime,namealt,namelat,namelon,))
bv_    = nc2.createVariable('n2','f', (nametime,namealt,namelat,namelon,))

###

alts_[:] = alts[:]
alts_.info = 'Altitude (m) m added by script createVenus_ex.py'
speed_[:] = speed[:]
speed_.info = 'Wind speed in m/s added by script createVenus_ex.py'
rot_[:] = rot[:]
rot_.info = 'Vorticity in s-1 added by script createVenus_ex.py'
#div_[:] = div[:]
#div_.info = 'Divergence in s-1 added by script createVenus_ex.py'
teta_[:] = teta[:]
teta_.info = 'Potential temperature added by script createVenus_ex.py'
#pv_[:] = pv[:]
#pv_.info = 'Ertel Potential Vorticity added by script addPV.py'
realw_[:] = realw[:]
realw_.info = 'Vertical velocity in m/s added by script createVenus_ex.py'
#bv_[:] = bv2[:]
#bv_.info = 'Square of Brunt Vaisala Frequency s-2 added by script createVenus_ex.py'


nc2.close()



