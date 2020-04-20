"""
Garrett Witzke
Constrain CAO statistics 
Hourly Contour Plots of LWP
XZ Slices
"""

import numpy as np
import struct
import netCDF4
from matplotlib import pylab as plt


with open('constrain.ini') as f:
    for line in f:
        if line.split('=')[0] == 'ktot':
            kmax = int(line.split('=')[1])
        if line.split('=')[0] == 'zsize':
            zsize = float(line.split('=')[1])

#----------------------------Constants----------------------------------------------
Rd = 287.0                #Universal Gas constant for dry air (J/K mol)
cp = 1006.0               #Specific Heat of Dry Air (J/kg K)
Lv = 2.26e6               #Latent Heat of Vaporization (J/kg)
p0 = 100900               #Surface pressure (Pa)
f  = 1.3e-4               #1/seconds Coriolis Parameter at a latitude of 63
Nc = 50                   #1/cm^3 Cloud Droplet conc

#----------------------------import results-----------------------------------------
stats  = netCDF4.Dataset('ql.nc', 'r')
stats1 = netCDF4.Dataset('qr.nc', 'r')

t_hourly = np.array([0,2,4,6,8,10,12,14,16,18,20,22,24,26,28])

#Contour Plot for LWP(ql_path) but for now its just ql
t = stats.variables["time"][:]
z = stats.variables["z"][:]
y = stats.variables["y"][:]
x = stats.variables["x"][:]
ql = stats.variables["ql"][:,:,:,:]

X,Z = np.meshgrid(x,z)
ql_path = np.zeros(t_hourly.size)

for t in range(t_hourly.size):
    ql_path = np.sum(ql[t,:,:,:],axis=2)*kmax
    fig, ax = plt.subplots()
    

    CS = ax.contourf(X,Z,ql_path, levels=100, cmap='gnuplot2')
    ax.set_xlabel("x-direction(m)")
    ax.set_ylabel("z-direction(m)")
    ax.set_title('Hour = ')

    fig.suptitle(r'$Liquid\/Water\/Path\/(kg/m^2)$')
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
    fig.colorbar(CS,cax=cbar_ax)
    fig.colorbar(CS, cax=cbar_ax)
        
plt.show()

    