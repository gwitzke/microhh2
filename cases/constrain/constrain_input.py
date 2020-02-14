"""
Garrett Witzke
CONSTRAIN CAO case setup using initial conditions
provided by the Met Office

outline
No ice microphysics...yet...
"""

import numpy as np 
import netCDF4 as nc 

float_type = "f8"

#Important Constants
Rd = 287.0
cp = 1005.0
Lv = 2.5e6
p0 = 1009
f  = 1.3e-4


# Get number of vertical levels and size form constrain.ini
with open('constrain.ini') as f:
    for line in f:
        if line.split('=')[0] == 'ktot':
            kmax = int(line.split('=')[1])
        if line.split('=')[0] == 'zsize':
            zsize = float(line.split('=')[1])

dz = zsize/kmax

#importing dimensions, time varying data, and initial conditions from met office
stats = netCDF4("constrain_setup_forcing.nc", "r")
t   = stats.variables["time"][:]        #Time values (s)
z   = stats.variables["z"][:]           #Height values (m)
sst = stats.variables["SST"][:]         #Time varying Sea Surface Temperature (K)
#lat = stats.variables["lat"][:,:]      #Height-Time Varying Latitude (degrees North)
#lon = stats.variables["lon"][:,:]      #Height -Time Varying Longitude (degrees East)
T   = stats.variables["T"][0,:]         #Initial Temperature as height increase (K)
thl = stats.variables["theta_l"][0,:]   #Initial Theta_l as height increases (K)
qt  = stats.variables["qt"][0,:]        #Initial qt as height increases (kg/kg)
ql  = stats.variables["qc"][0,:]        #Initial ql as height increases (kg/kg)
qv  = stats.variables["qv"][0,:]        #Initial qv as height increases (kg/kg)
u   = stats.variables["u"][0,:]         #Initial u-component of velocity as height increase (m/s)
v   = stats.variables["v"][0,:]         #Initial v-component of velocity as height increases (m/s)
vg = -15.0 - 0.0024*z                   #Initial Vgeo as height increases   (m/s)
ug = 0*z                                #Initial Ugeo as height increases   (m/s)

#Convert Units to SI



time_l = t[0:30:2]                     #Look at time intervals of every hour
time_ls = np.append(time_l, 52200)     #Add-in the last half hour to list


#Save all the input data to netCDF file
nc_file = nc.Dataset("constrain_input".nc",mode="w", datamodel="NETCDF4",clobber=False)

#Creat a Dimension for the height
nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))
nc_z[:] = z[:]

#Create a group called init for the initial profiles
nc_group_init = nc_file.createGroup("init")

nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt  = nc_group_init.createVariable("qt" , float_type, ("z"))
nc_u   = nc_group_init.createVariable("u"  , float_type, ("z"))
nc_ug  = nc_group_init.createVariable("ug" , float_type, ("z"))
nc_v   = nc_group_init.createVariable("v"  , float_type, ("z"))
nc_vg  = nc_group_init.createVariable("vg" , float_type, ("z"))

nc_thl[:] = thl[:]
nc_qt [:] = qt [:]
nc_u  [:] = u  [:]
nc_ug [:] = ug [:]
nc_v  [:] = v  [:]
nc_vg [:] = vg [:]

#Create a group call "timedep" for the time dependant variables
nc_group_timedep = nc_file.createGroup("timedep")

nc_group_timedep.createDimension("time", t.size)
nc_group_timedep.createDimension("time_ls", t.size)


nc_time      = nc_group_timedep.createVariable("time", time.size)
nc_thl_sbot  = nc_group_timedep.createVariable() #thl profile as affected by the SST
nc_qt_sbot   = nc_group_timedep.createVariable() #qt profile as affected by the SST
nc_sst       = nc_group_timedep.createVariable("sst" , float_type, ("time"))
#nc_lat = nc_group_timedep.createVariable("lat" , float_type, ("time)"))
#nc_lon = nc_group_timedep.createVariable("lon" , float_type, ("time"))

nc_time_surface[:] = time_surface[:]
nc_thl_sbot    [:] = sbotthl     [:]
nc_qt_sbot     [:] = sbotqt      [:]
nc_sst         [:] = sst         [:]


nc_file.close()
