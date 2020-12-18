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

#Importing initial conditions data from Met Office file
stats = netCDF4.Dataset("constrain_setup_forcing.nc", "r")


t       = stats.variables["time"][:]              #time
z       = stats.variables["z"][:]                 #Altitude
sst     = stats.variables["SST"][:]               #Sea Surface Temperature
wsubs   = stats.variables["wsubs"][:,:]           #Large Scale Vertical Velocity
u       = stats.variables["U"][:,:]               #East-West Wind Velocity
v       = stats.variables["V"][:,:]               #North-South Wind Velocity
T       = stats.variables["T"][:,:]               #Absolute Temperature
theta_l = stats.variables["theta_l"][:,:]         #Liquid Water Potential Temperature
qt      = stats.variables["qt"][:,:]              #Total Water Mixing Ratio
ql      = stats.variables["qc"][:,:]              #Liquid Water Mixing Ratio
Vgeo    = -15-0.0024*z
Ugeo    = 0*z
thr = t/3600            #Converting time into hours




#Plotting Time-Varying Sea-Surface Temperature and Subsidence at various heights
fig, axs = plt.subplots(2,1, 
                figsize = (8,8),
                sharex = True,
                tight_layout = True)

#Sea Surface Temperature
axs[0].plot(thr,sst)
axs[0].set_ylabel("SST [K]")
axs[0].set_yticks([276,278,280,282,284])
axs[0].set_ylim([276,284])
axs[0].text(0.5,283,'(a)')

#Large Scale Subsidence
axs[1].plot(thr,wsubs[:,12], '-', label = "z = 513 m", color = 'red')
axs[1].plot(thr,wsubs[:,21], ':', label = "z = 1533 m", color = 'blue')
axs[1].plot(thr,wsubs[:,27], '--',label = "z = 2513 m", color = 'black')
axs[1].set_xlabel("Time (hours)")
axs[1].set_ylabel("Subsidence [m/s]")
axs[1].set_xticks([0,5,10,15])
axs[1].set_xlim([0,15])
axs[1].set_yticks([-0.05,-0.04,-0.03, -0.02, -0.01, 0.00])
axs[1].set_ylim([-0.05,0.00])
axs[1].text(0.5,-0.005, '(b)')
axs[1].legend(loc='lower left')



fig, axs = plt.subplots(2,2,
                figsize = (8,8),
                tight_layout = True,
                sharey = True)

#Liquid Water Potential Temperature
axs[0,0].plot(theta_l[0,:32],z[:32]/1000, color = 'blue')           
axs[0,0].set_xlabel(r'$\theta_l\;\;[K]$')
axs[0,0].set_xticks([270,275,280,285])
axs[0,0].set_xlim([270,285])
axs[0,0].set_ylabel(r'$Height\;\;[km]$')
axs[0,0].set_yticks([0,1,2,3])
axs[0,0].set_ylim([0,3])
axs[0,0].text(270.5,2.8,'(a)')

#Water Content
axs[0,1].plot(qt[0,:32]*1000, z[:32]/1000, label = 'qt', color = 'blue')
axs[0,1].plot(ql[0,:32]*1000, z[:32]/1000, label = 'ql', color = 'red')
axs[0,1].set_xlabel(r'$[g/kg]$')
axs[0,1].set_xticks([0, 0.5, 1.0, 1.5, 2.0, 2.5])
axs[0,1].set_xlim([0,2.5])
axs[0,1].set_yticks([0,1,2,3])
axs[0,1].set_ylim([0,3])
axs[0,1].text(0.05, 2.8, '(b)')
axs[0,1].legend() 

#Zonal Winds
axs[1,0].plot(u[0,:32], z[:32]/1000, label ='U', color ='blue', ls = 'solid')
axs[1,0].plot(Ugeo[:32], z[:32]/1000, label = r'$U_{geo}$', color = 'red', ls = 'dashed')
axs[1,0].set_xlabel('[m/s]')
axs[1,0].set_xticks([-6,-4,-2,0,2,4,6])
axs[1,0].set_xlim([-6,6])
axs[1,0].set_ylabel(r'$Height\;\;[km]$')
axs[1,0].set_yticks([0,1,2,3])
axs[1,0].set_ylim([0,3])
axs[1,0].text(-5.5,2.8, '(c)')
axs[1,0].legend()

#Meridional Winds
axs[1,1].plot(v[0,:32], z[:32]/1000, label = 'V', color = 'blue', ls = 'solid')
axs[1,1].plot(Vgeo[:32], z[:32]/1000, label = r'$V_{geo}$', color = 'red', ls = 'dashed')
axs[1,1].set_xlabel('[m/s]')
axs[1,1].set_xticks([-25,-20,-15,-10])
axs[1,1].set_xlim([-25,-10])
axs[1,1].set_yticks([0,1,2,3])
axs[1,1].set_ylim([0,3])
axs[1,1].text(-24.5,2.8, '(d)')
axs[1,1].legend()

plt.show()