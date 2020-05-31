"""
Garrett Witzke
Constrain CAO statistics 
Contour Plots
"""

import numpy as np
import struct
import netCDF4
import math
from matplotlib import pylab as plt
from matplotlib import colors as colors

#----------------------------Constants----------------------------------------------
Rd = 287.0                #Universal Gas constant for dry air (J/K mol)
cp = 1006.0               #Specific Heat of Dry Air (J/kg K)
Lv = 2.26e6               #Latent Heat of Vaporization (J/kg)
p0 = 100900               #Surface pressure (Pa)
f  = 1.3e-4               #1/seconds Coriolis Parameter at a latitude of 63
Nc = 50                   #1/cm^3 Cloud Droplet conc


#----------------------------import results-----------------------------------------
default = netCDF4.Dataset("constrain.default.0000000.nc", 'r') 


#-------------------------------Contour Plots---------------------------------------------
#                               Figures 8 & 9
#                              De Roode et Al.
#                        t[6]= hour 3, t[24]= hour 12

stats1  = netCDF4.Dataset('w.nc', 'r')
stats2  = netCDF4.Dataset('u.nc', 'r')
stats3  = netCDF4.Dataset('v.nc', 'r')
stats4  = netCDF4.Dataset('T.nc', 'r')
stats5  = netCDF4.Dataset('ql.nc', 'r')
stats6  = netCDF4.Dataset('qr.nc', 'r')
stats7  = netCDF4.Dataset('thl.nc', 'r')

#Contourf plot for u
t = stats2.variables["time"][:]
z = stats2.variables["z"][:]
y = stats2.variables["y"][:]
xh = stats2.variables["xh"][:]
u = stats2.variables["u"][6,4,:,:]
u2 = stats2.variables["u"][24,4,:,:]

X2,Y2 = np.meshgrid(xh,y)

fig, (ax1,ax2) = plt.subplots(ncols=2, squeeze=True, sharey=True, figsize=(8,4))
fig.suptitle("u (m/s) @ 100m")

CS1 = ax1.contourf(X2,Y2,u,levels=100,cmap='gnuplot2')
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
CS2 = ax2.contourf(X2,Y2,u2,levels=100, cmap='gnuplot2')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
fig.colorbar(CS1,cax=cbar_ax)
fig.colorbar(CS2, cax=cbar_ax)
plt.savefig('u-contour.png')

#Contourf plot for v

t = stats3.variables["time"][:]
z = stats3.variables["z"][:]
yh = stats3.variables["yh"][:]
x = stats3.variables["x"][:]
v = stats3.variables["v"][6,4,:,:]
v2 = stats3.variables["v"][24,4,:,:]

X3,Y3 = np.meshgrid(x,yh)

fig, (ax1,ax2) = plt.subplots(ncols=2, squeeze=True, sharey=True, figsize=(8,4))
fig.suptitle("v (m/s) @ 100m")

CS1 = ax1.contourf(X3,Y3,v,levels=100,cmap='gnuplot2',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
CS2 = ax2.contourf(X3,Y3,v2,levels=100, cmap='gnuplot2')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
fig.colorbar(CS1,cax=cbar_ax)
fig.colorbar(CS2, cax=cbar_ax)
plt.savefig('v-contour.png')


#Contourf plot for w
t = stats1.variables["time"][:]
zh = stats1.variables["zh"][:]
y = stats1.variables["y"][:]
x = stats1.variables["x"][:]
w = stats1.variables["w"][6,4,:,:]
w2 = stats1.variables["w"][24,4,:,:]


X1,Y1 = np.meshgrid(x,y)

fig, (ax1,ax2) = plt.subplots(ncols=2, squeeze=True, sharey=True, figsize=(8,4))
fig.suptitle("w (m/s) @ 100m")

CS1 = ax1.contourf(X1,Y1,w,levels=100,cmap='gnuplot2',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
CS2 = ax2.contourf(X1,Y1,w2,levels=100, cmap='gnuplot2')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
fig.colorbar(CS1,cax=cbar_ax)
fig.colorbar(CS2, cax=cbar_ax)
plt.savefig('w-contour.png')
#Contour Plot of T

t = stats4.variables["time"][:]
z = stats4.variables["z"][:]
y = stats4.variables["y"][:]
x = stats4.variables["x"][:]
T = stats4.variables["T"][6,4,:,:]
T2 = stats4.variables["T"][24,4,:,:]

X4,Y4 = np.meshgrid(x,y)

fig, (ax1,ax2) = plt.subplots(ncols=2, squeeze=True, sharey=True, figsize=(8,4))
fig.suptitle("T (K) @ 100m")

CS1 = ax1.contourf(X4,Y4,T,levels=100,cmap='gnuplot2',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
CS2 = ax2.contourf(X4,Y4,T2,levels=100, cmap='gnuplot2')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
fig.colorbar(CS1,cax=cbar_ax)
fig.colorbar(CS2, cax=cbar_ax)
plt.savefig('T-contour.png')
#Contour Plot for qr

t = stats6.variables["time"][:]
z = stats6.variables["z"][:]
y = stats6.variables["y"][:]
x = stats6.variables["x"][:]
qr = stats6.variables["qr"][6,4,:,:]
qr2 = stats6.variables["qr"][24,4,:,:]

X6,Y6 = np.meshgrid(x,y)

fig, (ax1,ax2) = plt.subplots(ncols=2, squeeze=True, sharey=True, figsize=(8,4))
fig.suptitle("qr @ 100m")

CS1 = ax1.contourf(X6,Y6,qr,levels=100,cmap='gnuplot2',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
CS2 = ax2.contourf(X6,Y6,qr2,levels=100, cmap='gnuplot2')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
fig.colorbar(CS1,cax=cbar_ax)
fig.colorbar(CS2, cax=cbar_ax)
plt.savefig('qr-contour.png')

#Contour Plot for LWP(ql_path) but for now its just ql
t = stats5.variables["time"][:]
z = stats5.variables["z"][:]
y = stats5.variables["y"][:]
x = stats5.variables["x"][:]
ql = np.sum(stats5.variables["ql"][6,:,:,:],axis=0)
ql2 = np.sum(stats5.variables["ql"][24,:,:,:],axis=0)

X5,Y5 = np.meshgrid(x,y)

fig, (ax1,ax2) = plt.subplots(ncols=2, squeeze=True, sharey=True, figsize=(8,4))
fig.suptitle("ql @ 100m")

CS1 = ax1.contourf(X5,Y5,ql,levels=100,cmap='gnuplot2')
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
CS2 = ax2.contourf(X5,Y5,ql2,levels=100, cmap='gnuplot2')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
fig.colorbar(CS1,cax=cbar_ax)
fig.colorbar(CS2, cax=cbar_ax)
plt.savefig('ql-contour.png')


#Contour plot for thl
t = stats7.variables["time"][:]
z = stats7.variables["z"][:]
y = stats7.variables["y"][:]
x = stats7.variables["x"][:]
thl = stats7.variables["thl"][6,4,:,:]
thl2 = stats7.variables["thl"][24,4,:,:]



X7,Y7 = np.meshgrid(x,y)

fig, (ax1,ax2) = plt.subplots(ncols=2, squeeze=True, sharey=True, figsize=(8,4))
fig.suptitle("Theta_l @ 100m")

CS1 = ax1.contourf(X7,Y7,thl, levels=100, cmap='gnuplot2')
ax1.set_title("Time = Hour 3")
ax1.set_xlabel('x-direction (m)')
CS2 = ax2.contourf(X7,Y7, thl2, levels=100, cmap='gnuplot2')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
fig.subplots_adjust(right = 0.8)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
fig.colorbar(CS1, cax=cbar_ax)
fig.colorbar(CS2, cax=cbar_ax)
plt.savefig('thl-contour.png')


plt.show()


