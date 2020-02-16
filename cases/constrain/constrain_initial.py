"""
Garrett Witzke
Initial Conditions for the CONSTRAIN CAO case
Data provided by Met Office 
"http://appconv.metoffice.com/cold_air_outbreak/constrain_case/crm_setup.html"
"""
#Importing python packages
import numpy as np
import numpy.ma
import netCDF4
import struct 
from matplotlib import pylab as plt 

#Defining any functions i may need to do calculations
def find_nearest(array, value):                 #Finds value in array closest to user specified value
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


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
Ugeo    = 0*z


#Plotting Time-Varying Sea-Surface Temperature and Subsidence at various heights

thr = t/3600            #Converting time into hours

fig, axs = plt.subplots(2,1, figsize = (12,6))
axs[0].plot(thr,sst)
axs[0].set_title("Sea Surface Temperature")
axs[0].set_xlabel("Time (hours)")
axs[0].set_ylabel("Temperature (K)")
axs[0].axis([0,14.5,276,284])
axs[1].plot(thr,wsubs[:,12], '-', label = "z = 513 m")
axs[1].plot(thr,wsubs[:,21], ':', label = "z = 1533 m")
axs[1].plot(thr,wsubs[:,27], '--',label = "z = 2513 m")
axs[1].set_title("Large Scale Subsidence")
axs[1].set_xlabel("Time (hours)")
axs[1].set_ylabel("Subsidence")
axs[1].axis([0,14.5,-0.05,0])
axs[1].legend(loc='lower left')
fig.tight_layout()
plt.show()

#Plotting Liquid Water Potential Temp

plt.figure()
plt.plot(theta_l[0,:],z[:])           
plt.axis([270,285,0,3000])
plt.title("Liquid Water Potential Temperature")
plt.xlabel(r'$\theta_l (K)$')
plt.ylabel(r'$Height (m)$')
plt.show()

#Plotting Water Mixing Ratios
plt.figure()
plt.plot(qt[0,:]*1000,z, label = r'$q_t$')
plt.plot(ql[0,:]*1000,z, '--',label = r'$q_l$')
plt.axis([0,2.5,0,3000])
plt.title("Mixing Ratios")
plt.xlabel(r'$(g/kg)$')
plt.ylabel(r'$Height (m)$')
plt.legend(loc='upper right')
plt.show()

#Plotting of Vertical and Horizontal Wind Speeds Initial conditions

fig, axs = plt.subplots(1,2)
axs[0].plot(v[0,:],z[:],label = r'$V$')
axs[0].plot(Vgeo,z[:], '--',label = r'$V_{geo}$')
axs[0].set_xlabel("m/s")
axs[0].set_ylabel("height (m)")
axs[0].set_title("Vertical Wind Speed")
axs[0].axis([-25,-10,0,3000])
axs[0].legend(loc='lower left')
axs[1].plot(u[0,:],z[:], label = r'$U$')
axs[1].plot(Ugeo,z[:], '--', label = r'$U_{geo}$')
axs[1].set_xlabel("m/s")
axs[1].set_ylabel("height (m)")
axs[1].set_title("Horizontal Wind Speed")
axs[1].axis([-6,6,0,3000])
axs[1].legend(loc='lower left')
fig.tight_layout()
plt.show()

# Plotting Temperature and Temperature Gradient in Celsius

T[0,:] -= 273.15 #Convert to celsius
z[:] /= 1000     #Convert to kilometres

gamma = (-1) * np.gradient(T[0,:])          #Initial Temperature Lapse Rate defined as -dT/dz
a = find_nearest(gamma,2.01)                #Tropopause located at height where lapse rate = 2 and continues to decrease from then on in height

fig, axs = plt.subplots(2,1)
axs[0].plot(T[0,:],z[:])
axs[0].set_ylabel('Height (m)')
axs[0].set_xlabel(r'$T\/(^\circ C)$')
axs[0].set_title("Temperature")
axs[1].plot(gamma, z[:])
axs[1].set_ylabel('Height (m)')
axs[1].set_xlabel(r'$\Gamma_T\/(^\circ C / km)$')
axs[1].set_title("Temperature Gradient")
#axs[1].axvline(a,color = 'black')
#axs[1].axhline(y=8.25081, color = 'black')
axs[1].axis([-3,3,0,15])
fig.tight_layout()
plt.show()
