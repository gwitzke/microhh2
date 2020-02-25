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
Rd = 287.0                #Universal Gas constant for dry air
cp = 1005.0               #Specific Heat of Dry Air
Lv = 2.5e6                #Latent Heat of Vaporization 
p0 = 1009                 #Surface pressure hPa
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

#importing dimensions, time varying data, and initial conditions from met office
stats  = netCDF4.Dataset("constrain_setup_forcing.nc", "r")

t        = stats.variables["time"][:]          #Time values (s)
z        = stats.variables["z"][:38]           #Height values (m) (up to 5000m)
sst      = stats.variables["SST"][:]           #Time varying Sea Surface Temperature (K)
thl_sbot = stats.variables["SHF"][:]           #Time Varying surface forcing, sensible heat flux
qt_sbot  = stats.variables["LHF"][:]           #Time Varying surface forcing, latent heat flux
lat      = stats.variables["lat"][:,0]         #Time Varying Latitude (degrees North)
long     = stats.variables["lon"][:,0]         #Time Varying Longitude (degrees East)
wls      = stats.variables["wsubs"][:,:38]     #Initial Large Scale Vertical Velocity (m/s)
thl      = stats.variables["theta_l"][0,:38]   #Initial Theta_l as height increases (K)
thl_ls   = stats.variables["theta_l"][:,:38] #Initial large scale theta_l
qt       = stats.variables["qt_adj"][0,:38]    #Initial qt as height increases (kg/kg)
qt_ls    = stats.variables["qt_adj"][:,:38]  #Initial large scale qt
ql       = stats.variables["qc"][0,:38]        #Initial ql as height increases (kg/kg)
qv       = stats.variables["qv"][0,:38]        #Initial qv as height increases (kg/kg)
u        = stats.variables["U"][0,:38]         #Initial u-component of velocity as height increase (m/s)
v        = stats.variables["V"][0,:38]         #Initial v-component of velocity as height increases (m/s)
v_geo    = -15.0 - 0.0024*z                    #Initial Vgeo as height increases   (m/s)
u_geo    = 0*z                                 #Initial Ugeo as height increases   (m/s)

z_ls = np.array([2.50,      13.33,      33.33,      60.00,      93.33,
   133.33,     180.00,     233.33,     293.33,     360.00,
   433.33,     513.33,     600.00,     693.33,     793.33,
   900.00,    1013.33,    1133.33,    1260.00,    1393.33,
  1533.33,    1680.00,    1833.33,    1993.33,    2160.00,
  2333.33,    2513.33,    2700.00,    2893.33,    3093.33,
  3300.00,    3513.33,    3733.33,    3960.00,    4193.33,
  4433.33,    4680.00,    4933.33])                             #Large Scale Forcing Heights


time_ls = np.array([ 0.0, 3600.0, 7200.0, 10800.0, 14400.0,         
  18000.0, 21600.0, 25200.0, 28800.0, 32400.0,
  36000.0, 39600.0, 43200.0, 46800.0, 50400.0])                 #Large Scale Forcing time


#Save all the input data to netCDF file
nc_file = netCDF4.Dataset("constrain_input.nc",mode="w", datamodel="NETCDF4",clobber=False)

#Create a Dimension for the height
nc_file.createDimension("z", kmax)
nc_z   = nc_file.createVariable("z"  , float_type, ("z"))
nc_z   [:] = z[:]


#Create a group called init for the initial profiles
nc_group_init = nc_file.createGroup("init")

nc_thl    = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt     = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_u      = nc_group_init.createVariable("u"     , float_type, ("z"))
nc_ug     = nc_group_init.createVariable("u_geo" , float_type, ("z"))
nc_v      = nc_group_init.createVariable("v"     , float_type, ("z"))
nc_vg     = nc_group_init.createVariable("v_geo" , float_type, ("z"))


nc_thl      [:] = thl   [:]
nc_qt       [:] = qt    [:]
nc_u        [:] = u     [:]
nc_ug       [:] = u_geo [:]
nc_v        [:] = v     [:]
nc_vg       [:] = v_geo [:]



#Create a Dimension for time and large scale time
nc_file.createDimension("time", t.size)
nc_file.createDimension("time_ls" , time_ls.size)

#Create a group call "timedep" for the time dependant variables
nc_group_timedep = nc_file.createGroup("timedep") 

nc_time          = nc_group_timedep.createVariable("time"     , float_type, ("time"))
nc_time_ls       = nc_group_timedep.createVariable("time_ls"  , float_type, ("time_ls"))
nc_time      [:] = t       [:]
nc_time_ls   [:] = time_ls [:]

nc_sst           = nc_group_timedep.createVariable("sst"      , float_type, ("time"))
nc_thl_sbot      = nc_group_timedep.createVariable("thl_sbot" , float_type, ("time"))
nc_qt_sbot       = nc_group_timedep.createVariable("qt_sbot"  , float_type, ("time"))
nc_lat           = nc_group_timedep.createVariable("lat"      , float_type, ("time"))
nc_long          = nc_group_timedep.createVariable("long"     , float_type, ("time"))

nc_sst       [:] = sst [:]
nc_thl_sbot  [:] = thl_sbot [:]
nc_qt_sbot   [:] = qt_sbot  [:]
nc_lat       [:] = lat [:]
nc_long      [:] = long[:]

#Large Scale time dependent variables 

nc_thl_ls   = nc_group_timedep.createVariable("thl_ls" , float_type, ("time_ls", "z"))
nc_qt_ls    = nc_group_timedep.createVariable("qt_ls"  , float_type, ("time_ls", "z"))
nc_w_ls     = nc_group_timedep.createVariable("w_ls"   , float_type, ("time_ls", "z"))

nc_thl_ls  [:,:] = thl_ls  [:15,:]
nc_qt_ls   [:,:] = qt_ls   [:15,:]
nc_w_ls    [:,:] = wls     [:15,:]


nc_file.close()
 