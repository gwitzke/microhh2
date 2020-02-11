"""
Garrett Witzke
Initial Conditions for the CONSTRAIN CAO case
Data provided by Met Office 
"http://appconv.metoffice.com/cold_air_outbreak/constrain_case/crm_setup.html"
"""

import numpy as np
import numpy.ma
import netCDF4
import struct 
from matplotlib import pylab as plt 

#Importing initial conditions data from Met Office file
stats = netCDF4.Dataset("constrain_setup_forcing.nc", "r")
#array(time,z)
#array(time)

t = stats.variables["time"][:]                #time
z = stats.variables["z"][:]                   #Altitude
lat = stats.variables["lat"][:,:]             #Latitude
lon = stats.variables["lon"][:,:]             #Longitude
sst = stats.variables["SST"][:]               #Sea Surface Temperature
LHF = stats.variables["LHF"][:]               #Surface Latent Heat Flux
SHF = stats.variables["SHF"][:]               #Surface Sensible Heat Flux
wsubs = stats.variables["wsubs"][:,:]         #Large Scale Vertical Velocity
p = stats.variables["pressure"][:,:]          #Atmospheric Pressure
u = stats.variables["U"][:,:]                 #East-West Wind Velocity
v = stats.variables["V"][:,:]                 #North-South Wind Velocity
T = stats.variables["T"][:,:]                 #Absolute Temperature
theta = stats.variables["theta"][:,:]         #Potential Temperature
theta_l = stats.variables["theta_l"][:,:]     #Liquid Water Potential Temperature
qt = stats.variables["qt"][:,:]               #Total Water Mixing Ratio
qt_adj = stats.variables["qt_adj"][:,:]       #Total Water Mixing Ratio adjusted so initial profile saturated in cloud region
qv = stats.variables["qv"][:,:]               #Water Vapor Mixing Ratio
qv_adj = stats.variables["qv_adj"][:,:]       #Water Vapor Mixing Ratio adjusted so initial profile saturated in cloud region
ql = stats.variables["qc"][:,:]               #Liquid Water Mixing Ratio
RHw = stats.variables["RHw"][:,:]             #Initial relative humidity wrt water, calculated using theta_l and qt
RHi = stats.variables["RHi"][:,:]             #Initial relative humidity wrt ice, calculated using theta_l and qt
o3 = stats.variables["o3"][:]                 #Ozone Mass Mixing Ratio
LW = stats.variables["LWhr"][:,:]             #Long Wave Heating Rates
SW = stats.variables["SWhr"][:,:]             #Short Wave Heating Rates
Vgeo = -15-0.0024*z
Ugeo = 0

#Plotting Initial Conditions
plt.figure(1) 
plt.plot(t/3600,sst)
plt.title("Sea Surface Temperature")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.show()


plt.figure(2)
plt.plot(theta_l[0,0:35],z[0:35])           
plt.title("Liquid Water Potential Temperature")
plt.xlabel(r'$\theta_l (K)$')
plt.ylabel("Height (m)")
plt.show()


plt.figure(3)
plt.plot(qt[0,:],z, label = r'$q_t$')
plt.plot(ql[0,:],z, '--',label = r'$q_l$')
plt.title("Mixing Ratios")
plt.xlabel(r'$(g/kg)$')
plt.ylabel(r'$Height (m)$')
plt.legend(loc='upper right')
plt.show()

plt.figure(4)
plt.plot(v[0:],z[:])
plt.title("Vertical Wind Speed")
plt.show(
    
"""
plt.figure(4)
plt.plot(t/3600,wsubs[:,1], '-',label="z = 513 m")
plt.plot(t/3600,wsubs[:,2], '.', label="z = 1533 m")
plt.plot(t/3600,wsubs[:,3], '--', label="z = 2513")
plt.title("Large Scale Subsidence")
plt.xlabel(r'$Time$')
plt.ylabel(r'$Subsidence (m/s)$')
plt.legend(loc='lower left')
plt.show()
"""
