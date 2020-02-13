"""
Garrett Witzke
CONSTRAIN CAO case setup using initial conditions
provided by the Met Office

outline
"""

import numpy as np 
import netCDF4 as nc 

float_type = "f8"

# Get number of vertical levels and size form constrain.ini
with open('constrain.ini') as f:
    for line in f:
        if line.split('=')[0] == 'ktot':
            kmax = int(line.split('=')[1])
        if line.split('=')[0] == 'zsize':
            zsize = float(line.split('=')[1])

dz = zsize/kmax

#Importing initial conditions from Met Office File
stats = nc.Dataset("constrain_setup_forcing.nc", 'r',)

z = stats.variables["z"][:]
thl = stats.variables["theta_l"][:,:]
qt = stats.variables["qt]"][:,:]

#Importing the time series and setting it to hours
time = stats.variables["time"][:]
time_surface = time/3600



#Unit Convertions 


#Save all the input data to netCDF file
nc_file = nc.Dataset("constrain_input".nc",mode="w", datamodel="NETCDF4",clobber=False)

#Creat a Dimension for the height
nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))
nc_z[:] = z[:]

#Create a group called init for the initial profiles
nc_group_init = nc_file.createGroup("init")

#Create a group call "timedep" for the time dependant variables
nc_group_timedep = nc_file.createGroup("timedep")

nc_file.close()
