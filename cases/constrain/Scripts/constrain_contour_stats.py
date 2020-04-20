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
default = netCDF4.Dataset("constrain_default_0000000.nc", 'r') 


#Only grabbing necessary data

#Assigning basic variable names
t    = default.variables["time"][:]                              #Time
z    = default.variables["z"][:]                                 #Height 
zh   = default.variables["zh"][:38]                              #Half Height
zi   = default.groups["thermo"].variables["zi"][:]               #Boundary layer depth
rho  = default.groups["thermo"].variables["rho"][:,:]            #Density
rhoh = default.groups["thermo"].variables["rhoh"][:,:]           #Half level Density
Tt   = default.groups["thermo"].variables["T"][:,:]              #Absolute Temperature

#Assigning fractional area contained in mask
areat  = default.groups["default"].variables["area"][:,:]
areaht = default.groups["default"].variables["areah"][:,:]

#Assigning Time Series Variable names
LWPt  = default.groups["thermo"].variables["ql_path"][:]*1000   #Liquid Water Path units of g/m^2
RWPt = default.groups["default"].variables["qr_path"][:]*1000   #Rain Water Specific Humidity Path g/m^2
nr_path  = default.groups["default"].variables["nr_path"][:]    #Number Density Rain Path 
rainflx = default.groups["thermo"].variables["rr"][:]           #Mean Surface Rain Rate
thlflux = default.groups["default"].variables["thl_flux"][:,:]     #turbulent flux of theta_l
qtflux = default.groups["default"].variables["qt_w"][:,:]       #Turbulent flux of qt                      
    
#Assigning Thermodynamic variable names
thl  = default.groups["default"].variables["thl"][:,:]        #Theta_l
qtt  = default.groups["default"].variables["qt"][:,:]*1000    #qt units of g/kg
qlt  = default.groups["thermo"].variables["ql"][:,:]*1000     #ql units of g/kg
qrt  = default.groups["default"].variables["qr"][:,:]*1000    #qr units of g/kg
ql_cover = default.groups["thermo"].variables["ql_cover"][:]
qvt  = qtt - qlt
ql_frac = default.groups["thermo"].variables["ql_frac"][:,:]


#Assigning Wind Component variable names
wt     = default.groups["default"].variables["w"][:,:]        #Vertical Velocity (z-direction)
ut     = default.groups["default"].variables["u"][:,:]        #U direction velocity (x- component)
u2t    = default.groups["default"].variables["u_2"][:,:]      
vt     = default.groups["default"].variables["v"][:,:]        #V direction velocity (y -component)
v2t    = default.groups["default"].variables["v_2"][:,:]
w2t    = default.groups["default"].variables["w_2"][:,:]
ufluxt = default.groups["default"].variables["u_w"][:,:]      #Turblulent flux of u velocity
vfluxt = default.groups["default"].variables["v_w"][:,:]      #Turbulent flux of v velocity
Ufluxt = ufluxt + vfluxt
tket   = 0.5*(u2t + v2t + 0.5*(w2t[:,0:-1] + w2t[:,1::]))     #Turbulent Kinetic Energy



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


#Contourf plot for u
t = stats2.variables["time"][:]
z = stats2.variables["z"][:]
y = stats2.variables["y"][:]
xh = stats2.variables["xh"][:]
u = stats2.variables["u"][6,5,:,:]
u2 = stats2.variables["u"][24,5,:,:]

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


#Contourf plot for v

t = stats3.variables["time"][:]
z = stats3.variables["z"][:]
yh = stats3.variables["yh"][:]
x = stats3.variables["x"][:]
v = stats3.variables["v"][6,5,:,:]
v2 = stats3.variables["v"][24,5,:,:]

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



#Contourf plot for w
t = stats1.variables["time"][:]
zh = stats1.variables["zh"][:]
y = stats1.variables["y"][:]
x = stats1.variables["x"][:]
w = stats1.variables["w"][6,5,:,:]
w2 = stats1.variables["w"][24,5,:,:]


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

#Contour Plot of T

t = stats4.variables["time"][:]
z = stats4.variables["z"][:]
y = stats4.variables["y"][:]
x = stats4.variables["x"][:]
T = stats4.variables["T"][6,5,:,:]
T2 = stats4.variables["T"][24,5,:,:]

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

#Contour Plot for qr

t = stats6.variables["time"][:]
z = stats6.variables["z"][:]
y = stats6.variables["y"][:]
x = stats6.variables["x"][:]
qr = stats6.variables["qr"][6,5,:,:]
qr2 = stats6.variables["qr"][24,5,:,:]

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

CS1 = ax1.contourf(X5,Y5,ql,levels=100,cmap='gnuplot2',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
CS2 = ax2.contourf(X5,Y5,ql2,levels=100, cmap='gnuplot2')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
fig.colorbar(CS1,cax=cbar_ax)
fig.colorbar(CS2, cax=cbar_ax)


plt.show()


