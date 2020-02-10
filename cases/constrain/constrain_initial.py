"""
Garrett Witzke
Initial Conditions for the CONSTRAIN CAO case
Data provided by Met Office 
"http://appconv.metoffice.com/cold_air_outbreak/constrain_case/crm_setup.html"
"""

import numpy
import numpy.ma
import netCDF4
import struct 
from matplotlib import pylab as plt 

#Importing initial conditions data from Met Office file
stats = netCDF4.Dataset("constrain_setup_forcing.nc", "r")

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
qt_adj = stats.variables["qt_adj"][:,:]       #Total Water Mixing Ratio adjusted so initial profeil saturated in cloud region
qv = stats.variables["qv"][:,:]               #Water Vapor Mixing Ratio
qv_adj = stats.variables["qv_adj"][:,:]       #Water Vapor Mixing Ratio adjusted so initial profile saturated in cloud region
ql = stats.variables["qc"][:,:]               #Liquid Water Mixing Ratio
RHw = stats.variables["RHw"][:,:]             #Initial relative humidity wrt water, calculated using theta_l and qt
RHi = stats.variables["RHi"][:,:]             #Initial relative humidity wrt ice, calculated using theta_l and qt
o3 = stats.variables["o3"][:]                 #Ozone Mass Mixing Ratio
LW = stats.variables["LWhr"][:,:]             #Long Wave Heating Rates
SW = stats.variables["SWhr"][:,:]             #Short Wave Heating Rates

#Testing Masked Array accessibility


#Plotting Initial Conditions
plt.figure(1) 
plt.plot(t,sst)
plt.title("Sea Surface Temperature")
plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")
plt.show()

plt.figure(2)
plt.plot()
plt.title()
plt.xlabel()
plt.ylabel()
plt.show()