"""
Garrett Witzke
CONSTRAIN CAO case setup using initial conditions
provided by the Met Office

outline
No ice microphysics...yet...
"""

import numpy as np 
import netCDF4 

float_type = "f8"

#Important Constants
Rd = 287.0                #Universal Gas constant for dry air (J/K mol)
cp = 1006.0               #Specific Heat of Dry Air (J/kg K)
Lv = 2.26e6                #Latent Heat of Vaporization (J/kg)
p0 = 100900               #Surface pressure (Pa)
f  = 1.3e-4               #1/seconds Coriolis Parameter at a latitude of 63
Nc = 50                   #1/cm^3 Cloud Droplet concentration number


# Get number of vertical levels and size from constrain.ini
with open('constrain.ini') as f:
    for line in f:
        if line.split('=')[0] == 'ktot':
            kmax = int(line.split('=')[1])
        if line.split('=')[0] == 'zsize':
            zsize = float(line.split('=')[1])

dz = zsize/kmax
height = 38 #Index slice to sets max z height to 5000m (roughly)

#importing dimensions, time varying data, and initial conditions from met office
stats  = netCDF4.Dataset("constrain_setup_forcing.nc", "r")

t        = stats.variables["time"][:]                   #Time values (s)
z        = stats.variables["z"][:height]                #Height values (m) (up to 5000m)
sst      = stats.variables["SST"][:]                    #Time varying Sea Surface Temperature (K)
SHF      = stats.variables["SHF"][:]                    #Time Varying surface forcing, sensible heat flux
LHF      = stats.variables["LHF"][:]                    #Time Varying surface forcing, latent heat flux
lat      = stats.variables["lat"][:,0]                  #Time Varying Latitude (degrees North)
lon      = stats.variables["lon"][:,0]                  #Time Varying Longitude (degrees East)
LWhr     = stats.variables["LWhr"][0,:height]           #Initial Longwave Radiation
SWhr     = stats.variables["SWhr"][0,:height]           #Initial Shortwave Radiation
w_ls     = stats.variables["wsubs"][:,:height]          #Initial Large Scale Vertical Velocity (m/s)
thl      = stats.variables["theta_l"][0,:height]        #Initial Theta_l as height increases (K)
qt       = stats.variables["qt_adj"][0,:height]         #Initial qt as height increases (kg/kg)  
qv       = stats.variables["qv_adj"][0,:height]         #Initial qv as height increases (kg/kg)
u        = stats.variables["U"][0,:height]              #Initial u-component of velocity as height increase (m/s)
v        = stats.variables["V"][0,:height]              #Initial v-component of velocity as height increases (m/s)
v_geo    = -15.0 - 0.0024*z                             #Initial Vgeo as height increases   (m/s)
u_geo    = 0*z                                          #Initial Ugeo as height increases   (m/s)

#Creating empty variables for qr and nr, whatever those are...
qr = np.zeros(z.size)
nr = np.zeros(z.size)

#Calculating the surface fluxes
rho = p0/(Rd*thl[0]*(1. + 0.61*qt[0]))
sbotthl = SHF/(rho*cp)
sbotqt  = LHF/(rho*Lv)

#Save all the input data to netCDF file
nc_file = netCDF4.Dataset("constrain_input.nc",mode="w", datamodel="NETCDF4",clobber=False)

#Initial Conditions

#Create a Dimension for the height
nc_file.createDimension("z", kmax)
nc_z   = nc_file.createVariable("z"  , float_type, ("z"))

nc_z   [:] = z[:]

#Create a group called init for the initial profiles
nc_group_init = nc_file.createGroup("init")

nc_thl    = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt     = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_qv     = nc_group_init.createVariable("qv"    , float_type, ("z"))
nc_u      = nc_group_init.createVariable("u"     , float_type, ("z"))
nc_ug     = nc_group_init.createVariable("u_geo" , float_type, ("z"))
nc_v      = nc_group_init.createVariable("v"     , float_type, ("z"))
nc_vg     = nc_group_init.createVariable("v_geo" , float_type, ("z"))
nc_qr     = nc_group_init.createVariable("qr"    , float_type, ("z"))
nc_nr     = nc_group_init.createVariable("nr"    , float_type, ("z"))


nc_thl      [:] = thl   [:]
nc_qt       [:] = qt    [:]
nc_qv       [:] = qv    [:]
nc_u        [:] = u     [:]
nc_ug       [:] = u_geo [:]
nc_v        [:] = v     [:]
nc_vg       [:] = v_geo [:]
nc_qr       [:] = qr    [:]
nc_nr       [:] = nr    [:]



#Create a Dimension for time and large scale time
nc_file.createDimension("time", t.size)

#Create a group call "timedep" for the time dependant variables
nc_group_timedep = nc_file.createGroup("timedep") 
nc_time          = nc_group_timedep.createVariable("time"     , float_type, ("time"))
nc_time      [:] = t       [:]


nc_sst           = nc_group_timedep.createVariable("sst"      , float_type, ("time"))
nc_thl_sbot      = nc_group_timedep.createVariable("thl_sbot" , float_type, ("time"))
nc_qt_sbot       = nc_group_timedep.createVariable("qt_sbot"  , float_type, ("time"))
nc_lat           = nc_group_timedep.createVariable("lat"      , float_type, ("time"))
nc_lon           = nc_group_timedep.createVariable("lon"      , float_type, ("time"))
nc_wls           = nc_group_timedep.createVariable("wls"      , float_type, ("time","z"))

nc_sst       [:] = sst     [:]
nc_thl_sbot  [:] = sbotthl [:]
nc_qt_sbot   [:] = sbotqt  [:]
nc_lat       [:] = lat     [:]
nc_lon       [:] = lon     [:]
nc_wls       [:,:] = w_ls   [:,:]





nc_file.close()
 
