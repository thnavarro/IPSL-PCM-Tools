#!/usr/bin/env python


from    optparse   import  OptionParser
import  numpy      as      np
from    netCDF4    import  Dataset
import  gc


# createVenus_airglow:
# create a file with O2 airglow from a Venus (LMD) GCM output.
# from Gabriella email 10 March 2020

######################
############ Constants
nu = 0.35
cp0 = 1000.
pref = 9.2e6
temp0 = 460. 
R = 191.4 # oui mais non car la composition change
R0 = 8.3144626
NAVO = 6.0221367E+23 # [mol-1]
RVenus = 6051.8   # [km]
Mm = 43.44
kboltzman = 1.381e-23  # [J/K]
g = 8.87
######################
kappa0 = R/cp0 
######################


## compute vertical vs horizontal quantities ?
# Preliminary
#tracer = True
#tracer = False

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

#############################
### min level
lmin = 38 # 1000 Pa, very safe
ap = ap[lmin:]
bp = bp[lmin:]
############################

parser = OptionParser()
parser.add_option('-t','--tracer', action='store_true',        dest='tracer',default=False, help="compute tracer quality of O2 singlet delta")
parser.add_option('-l','--level',  action='store', type='int', dest='level', default=1,     help="output level option [1: 2D variables] [2: 3D variables] [3: more 3D variables for debug]")
parser.add_option('--co',          action='store_true',        dest='co',    default=False, help="compute CO")


(opt,args) = parser.parse_args()
print(opt.level)

if len(args)==0:
    print('Stop here! I need a file as argument(s)!'); exit()
else:
    fname = args[0]

nc = Dataset(fname,'r+')
ncv = nc.variables


lat  = ncv['lat'][:]
lon  = ncv['lon'][:]
alt  = ncv['presnivs'][lmin:]
time  = ncv['time_counter'][:]

#it = np.arange(60)
it = np.arange(len(time))
time = time[it]

o_vmr    = ncv['o'][it,lmin:,:,:]
#o2_vmr    = ncv['o2'][it,lmin:,:,:]
co2_vmr    = ncv['co2'][it,lmin:,:,:]
temp     = ncv['temp'][it,lmin:,:,:]

if opt.tracer:
  vitu = ncv['vitu'][it,lmin:,:,:]
  vitv = ncv['vitv'][it,lmin:,:,:]
  vitw = ncv['vitw'][it,lmin:,:,:]
  mmean= ncv['mmean'][it,lmin:,:,:]
  o2dg_vmr = ncv['o2dg'][it,lmin:,:,:]

co = False
if 'co' in ncv and opt.co:
  co = True
  co_vmr = ncv['co'][it,lmin:,:,:]
  co_vmr[:,alt>25.,:,:] = 0.

aire = ncv['aire'][:,:]
tops = ncv['tops'][it,:,:]



#########################################
######### Compute density

#if 'pres' in ncv:
if 'plop' in ncv:
  pres = ncv['pres'][it,:,:,:]
else:
  ps = ncv['psol'][it,:,:]

  pstiled = ps[:,np.newaxis,:,:]
  aptiled = ap[np.newaxis,:,np.newaxis,np.newaxis]
  bptiled = bp[np.newaxis,:,np.newaxis,np.newaxis]

  pres_interlay = pstiled*bptiled + aptiled
  gc.collect()
  #pres = 0.5*(pres_interlay[:,1:len(ap),:,:] + pres_interlay[:,0:len(ap)-1,:,:])
  pres = 0.5*(pres_interlay[:,1:,:,:] + pres_interlay[:,:-1,:,:])
  alts = -R*temp*np.diff(np.log(pres_interlay),axis=1)/g # level thickness in meters
  del pres_interlay,pstiled,aptiled,bptiled
  gc.collect()
  print('createVenus_airglow.py: 3D pressure reconstructed from surface pressure.')

alts[:,-1,:,:] = alts[:,-2,:,:]
print('Average top GCM altitude is '+str(np.sum(alts,axis=1).mean(0).mean(0).mean(0)*1e-3)+' km')

dens = pres*NAVO*1e-6/(R0*temp)       #  [molecule cm-3]

cc_o   = dens * o_vmr

#cc_o2dg  = dens * o2dg_vmr
#cc_co2  = dens * co2_vmr

del o_vmr
gc.collect()

#######################################
##
### Coefficient rates FROM F.Lefevre
##
#######################################

#a002 = 2.5*5.2e-35*np.exp(900./temp)*dens         # R1:  k from Tsang and hampson, 1986
a002 = 2.5*9.46e-34*np.exp(485./temp)*dens         # R1:  k from Cambpell and Grey 1973
epsilon = 0.75
tau  = 4470                                     # R2: 1/A  (A =einstein coeff.) s-1 =
quenching = 0.5e-20                             # R2: C   (value given by FL)

#####################################################################################################
##
### Coefficient rates from Gerard et al. 2014   (to be eventually used later for comparison with obs..)
##
######################################################################################################

#tau = 4566.21
#quenching = 2e-20
#quenching = 1.e-20                         # R2: Values used in Soret 2011



##################################################################
#    O2dg local Volume Emission Rate (VER)
#
#  USEFUL CONVERSIONS
#  10**4 photons.cm-3.s-1 = 1kR.km-1  (for the limb integration)
#
#  10**10 ph m-2 s-1       = 1R       (column emission)
#  10**12 ph cm-2 s-1      = 1MR
##################################################################


emiss_o2 = (a002*epsilon*cc_o*cc_o)/(1. + tau*quenching*dens)    # photons.cm-3.s-1
emiss_o2[emiss_o2>1e8]  = 0. # some weird high values with r2191

col_emiss_o2 = np.sum(emiss_o2*alts*1e2,axis=1) # [ph cm-2 s-1]

if opt.level < 2:
  del emiss_o2
  gc.collect()

col_o = cc_o*alts*1e2 # [molecule cm-2]
col_o[:,:70-lmin,:,:] = 0.
col_o = np.sum(col_o,axis=1) # [molecule cm-2]


if co:
  col_co = np.sum(dens*co_vmr*alts*1e2,axis=1) # [molecule cm-2]
  del co_vmr
  gc.collect()

tot_col_emiss_o2 = np.sum(np.sum(col_emiss_o2*aire[np.newaxis,:,:],axis=1),axis=1)/np.sum(np.sum(aire,axis=0),axis=0) # [ph cm-2 s-1]

day_col_o = col_o + 0.
day_col_o[tops!=0] = 0.
tot_day_col_o = np.sum(np.sum(day_col_o*aire[np.newaxis,:,:],axis=1),axis=1)/2 # molecule



# ############# Using GCM's O2dg
# emiss_o2_bis = (dens*o2dg_vmr)/(1. + tau*quenching*dens)    # photons.cm-3.s-1
# col_emiss_o2_bis = np.sum(emiss_o2_bis*alts*1e2,axis=1) # [ph cm-2 s-1]


if opt.tracer:

   # moles of singlet oxygen
   cc_odg   = dens*1e6*o2dg_vmr*alts*aire[np.newaxis,np.newaxis,:,:]/NAVO # [moles]
   # horizontal gradients
   dxq = np.diff(np.concatenate((cc_odg,cc_odg[:,:,:,0][:,:,:,np.newaxis]),axis=3),axis=3)
   dyq = np.diff(np.concatenate((cc_odg,cc_odg[:,:,0,:][:,:,np.newaxis,:]),axis=2),axis=2)
   # horizontal advection terms
   h1 = vitu*dxq*len(lon)/(RVenus*1e3*np.pi*2*np.cos(lat*np.pi/180)[np.newaxis,np.newaxis,:,np.newaxis])
   h2 = vitv*dyq*len(lat)/(RVenus*1e3*np.pi)
   
   h1[:,:,0,:] = 0.
   h1[:,:,-1,:] = 0.
   h2[:,:,0,:] = 0.
   h2[:,:,-1,:] = 0.
   
   #h = np.abs(np.sum(h1,1)) + np.abs(np.sum(h2,1))
   h3d = h1 + h2
   
   del h1,h2,vitu,vitv,dens,o2dg_vmr
   gc.collect()
   

   # vertical gradient
   dzq = -np.diff(np.concatenate((cc_odg,cc_odg[:,0,:,:][:,np.newaxis,:,:]),axis=1),axis=1)
   # vertical advection term
   v3d = vitw*dzq*R0*temp/(alts*mmean*1.e-3*g*pres)
   if opt.level == 3:
      realw = vitw*R0*temp/(mmean*1.e-3*g*pres)
   
   del vitw,mmean,temp,pres
   gc.collect()
   if opt.level != 3:
     del dzq,dxq,dyq 
     gc.collect()

##### d/dt chem

   # cell volume
   vol    = aire[np.newaxis,np.newaxis,:,:]*alts*1e6 # volume [cm3]
   # sink of atomic oxygen  = source of singlet oxygen
   source = a002*cc_o*cc_o*vol/NAVO  # cm-3.s-1.cm3/mol-1 = mol.s-1
   # sink of singlet oxygen with a 75 min time constant
   sink   = -2.2e-4*cc_odg # mol.s-1
##########################################
#### correction after review:
   source = epsilon*source
   # total moles is dens*vol/NAVO:
   sink   = sink - quenching*cc_odg*vol/NAVO
##########################################
##########################################
   c3d    = source + sink

   del alts,vol,a002,cc_o
   gc.collect()
   if opt.level < 3:
      del source,sink
      gc.collect()

   #ratio v/(h+v+c)
   x3d = np.abs(v3d)/(np.abs(h3d)+np.abs(v3d)+np.abs(c3d))
   #x   = np.mean(x3d,1)

   #ratio c/(h+v+c)
   y3d = np.abs(c3d)/(np.abs(h3d)+np.abs(v3d)+np.abs(c3d))
   #y   = np.mean(y3d,1)

   
else:

   del dens,temp,pres,a002,cc_o,alts
   gc.collect()


fnameout=fname[0:len(fname)-3]+'_airglow3.nc'

nc2 = Dataset(fnameout, 'w', format='NETCDF3_64BIT')
nc2.createDimension('lon',len(lon))
nc2.createDimension('lat',len(lat))
nc2.createDimension('presnivs',len(alt))
nc2.createDimension('time_counter',len(time))

lon_ = nc2.createVariable('lon', 'f', ('lon',))
lon_[:] = lon[:]
lat_ = nc2.createVariable('lat', 'f', ('lat',))
lat_[:] = lat[:]
alt_ = nc2.createVariable('presnivs', 'f', ('presnivs',))
alt_[:] = alt[:]
time_ = nc2.createVariable('time_counter', 'f', ('time_counter',))
time_[:] = time[:]


#alts_ = nc2.createVariable('alts','f', ('time_counter','presnivs','lat','lon',))
if opt.level >= 2:
  ag_ = nc2.createVariable('o2_emis','f', ('time_counter','presnivs','lat','lon',))
agcol_ = nc2.createVariable('o2_col','f', ('time_counter','lat','lon',))
#agcol2_ = nc2.createVariable('o2_col2','f', ('time_counter','lat','lon',))
ocol_ = nc2.createVariable('o_col','f', ('time_counter','lat','lon',))
if co:
  cocol_ = nc2.createVariable('co_col','f', ('time_counter','lat','lon',))
o2tot_ = nc2.createVariable('o2_tot','f', ('time_counter',))
otot_ = nc2.createVariable('o_tot','f', ('time_counter',))
if opt.tracer:
  #x_   = nc2.createVariable('x','f', ('time_counter','lat','lon',))
  x3d_ = nc2.createVariable('x3d','f', ('time_counter','presnivs','lat','lon',))
  #y_   = nc2.createVariable('y','f', ('time_counter','lat','lon',))
  y3d_ = nc2.createVariable('y3d','f', ('time_counter','presnivs','lat','lon',))
  if opt.level >= 3:
    h3d_ = nc2.createVariable('h3d','f', ('time_counter','presnivs','lat','lon',))
    v3d_ = nc2.createVariable('v3d','f', ('time_counter','presnivs','lat','lon',))
    c3d_ = nc2.createVariable('c3d','f', ('time_counter','presnivs','lat','lon',))
    source_ = nc2.createVariable('source','f', ('time_counter','presnivs','lat','lon',))
    sink_ = nc2.createVariable('sink','f', ('time_counter','presnivs','lat','lon',))
    dxq_ = nc2.createVariable('dxq','f', ('time_counter','presnivs','lat','lon',))
    dyq_ = nc2.createVariable('dyq','f', ('time_counter','presnivs','lat','lon',))
    dzq_ = nc2.createVariable('dzq','f', ('time_counter','presnivs','lat','lon',))
    realw_ = nc2.createVariable('realw','f', ('time_counter','presnivs','lat','lon',))
###
#alts_[:] = alts[:]
if opt.level >= 2:
  ag_[:] = emiss_o2[:]*1e-4
  ag_.info = 'O2 airglow VER in kR.km-1'
agcol_[:] = col_emiss_o2[:]*1e-12
agcol_.info = 'O2 airglow column emission in MR'
#agcol2_[:] = col_emiss_o2_bis[:]*1e-12
#agcol2_.info = 'O2 airglow column emission in MR'
ocol_[:] = col_o[:]
ocol_.info = 'O column in molecule cm-2'
if co:
  cocol_[:] = col_co[:]
  cocol_.info = 'CO column in molecule cm-2'
o2tot_[:] = tot_col_emiss_o2[:]*1e-12
o2tot_.info = 'Total O2 airglow column emission in MR'
otot_[:] = tot_day_col_o[:]
otot_.info = 'Total dayside O molecules'
if opt.tracer:
  #x_[:]   = x[:]
  #x_.info = 'v/(h+v+c)'
  x3d_[:] = x3d[:]
  x3d_.info = 'v/(h+v+c)'
  #y_[:]   = y[:]
  #y_.info = 'c/(h+v+c)'
  y3d_[:] = y3d[:]
  y3d_.info = 'c/(h+v+c)'
  if opt.level >= 3:
    h3d_[:] = h3d[:]
    v3d_[:] = v3d[:]
    c3d_[:] = c3d[:]
    source_[:] = source[:]
    sink_[:] = sink[:]
    dxq_[:] = dxq[:]
    dyq_[:] = dyq[:]
    dzq_[:] = dzq[:]
    realw_[:] = realw[:]


nc2.close()



