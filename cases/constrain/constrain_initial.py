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


t       = stats.variables["time"][:]              #time
z       = stats.variables["z"][:]                 #Altitude
lat     = stats.variables["lat"][:,:]             #Latitude
lon     = stats.variables["lon"][:,:]             #Longitude
sst     = stats.variables["SST"][:]               #Sea Surface Temperature
LHF     = stats.variables["LHF"][:]               #Surface Latent Heat Flux
SHF     = stats.variables["SHF"][:]               #Surface Sensible Heat Flux
wsubs   = stats.variables["wsubs"][:,:]           #Large Scale Vertical Velocity
p       = stats.variables["pressure"][:,:]        #Atmospheric Pressure
u       = stats.variables["U"][:,:]               #East-West Wind Velocity
v       = stats.variables["V"][:,:]               #North-South Wind Velocity
T       = stats.variables["T"][:,:]               #Absolute Temperature
theta   = stats.variables["theta"][:,:]           #Potential Temperature
theta_l = stats.variables["theta_l"][:,:]         #Liquid Water Potential Temperature
qt      = stats.variables["qt"][:,:]              #Total Water Mixing Ratio
qt_adj  = stats.variables["qt_adj"][:,:]          #Total Water Mixing Ratio adjusted so initial profile saturated in cloud region
qv      = stats.variables["qv"][:,:]              #Water Vapor Mixing Ratio
qv_adj  = stats.variables["qv_adj"][:,:]          #Water Vapor Mixing Ratio adjusted so initial profile saturated in cloud region
ql      = stats.variables["qc"][:,:]              #Liquid Water Mixing Ratio
RHw     = stats.variables["RHw"][:,:]             #Initial relative humidity wrt water, calculated using theta_l and qt
RHi     = stats.variables["RHi"][:,:]             #Initial relative humidity wrt ice, calculated using theta_l and qt
o3      = stats.variables["o3"][:]                #Ozone Mass Mixing Ratio
LW      = stats.variables["LWhr"][:,:]            #Long Wave Heating Rates
SW      = stats.variables["SWhr"][:,:]            #Short Wave Heating Rates
Vgeo    = -15-0.0024*z
Ugeo    = 0 *z

#Converting Units
thr = t/3600            #Time in hours

plt.figure(1,figsize = (12,6)) 
plt.plot(thr,sst)
plt.axis([0,15,276,284])
plt.title("Sea Surface Temperature")
plt.xlabel(r'$Time (s)$')
plt.ylabel(r'$Temperature (K)$')
plt.show()


plt.figure(2, figsize = (12,6))
plt.plot(thr,wsubs[:,12], '-',label="z = 513 m")            
plt.plot(thr,wsubs[:,21], ':', label="z = 1533 m")          
plt.plot(thr,wsubs[:,27], '--', label="z = 2513 m")         
plt.axis([0,15,-0.05,0])
plt.title("Large Scale Subsidence")
plt.xlabel(r'$Time$')
plt.ylabel(r'$Subsidence (m/s)$')
plt.legend(loc='lower left')
plt.show()


plt.figure(3)
plt.plot(theta_l[0,:],z[:])           
plt.axis([270,285,0,3000])
plt.title("Liquid Water Potential Temperature")
plt.xlabel(r'$\theta_l (K)$')
plt.ylabel(r'$Height (m)$')
plt.show()


plt.figure(4)
plt.plot(qt[0,:]*1000,z, label = r'$q_t$')
plt.plot(ql[0,:]*1000,z, '--',label = r'$q_l$')
plt.axis([0,2.5,0,3000])
plt.title("Mixing Ratios")
plt.xlabel(r'$(g/kg)$')
plt.ylabel(r'$Height (m)$')
plt.legend(loc='upper right')
plt.show()


plt.figure(5)
plt.plot(v[0,:],z[:], label = r'$V$')
plt.plot(Vgeo,z, '--', label = r'$V_{geo}$')
plt.axis([-25,-10,0,3000])
plt.title("Vertical Wind Speed")
plt.xlabel(r'$(m/s)$')
plt.ylabel(r'$Height(m)$')
plt.legend(loc="upper right")
plt.show()


plt.figure(6)
plt.plot(u[0,:],z[:], label = r'$U$')
plt.plot(Ugeo,z, '--', label = r'$U_{geo}$')
plt.axis([-6,6,0,3000])
plt.title("Horizontal Wind Speed")
plt.xlabel(r'$(m/s)$')
plt.ylabel(r'$height(m)$')
plt.legend(loc='upper left')
plt.show()


plt.figure(7)
plt.plot(T[0,:], z[:])
plt.title("Temperature")
plt.xlabel(r'$Temperature (K)$')
plt.ylabel(r'$height (m)$')
plt.show()
