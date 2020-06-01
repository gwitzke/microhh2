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
stats1  = netCDF4.Dataset('w.nc', 'r')
stats2  = netCDF4.Dataset('u.nc', 'r')
stats3  = netCDF4.Dataset('v.nc', 'r')
stats4  = netCDF4.Dataset('T.nc', 'r')
stats5  = netCDF4.Dataset('ql.nc', 'r')
stats6  = netCDF4.Dataset('qr.nc', 'r')
stats7  = netCDF4.Dataset('thl.nc', 'r')

#-------------------------------Contour Plots---------------------------------------------
#                               Figures 8 & 9
#                              De Roode et Al.
#                        t[6]= hour 3, t[24]= hour 12



#Liquid Water Path

t = stats5.variables["time"][:]
z = stats5.variables["z"][:]
y = stats5.variables["y"][:]
x = stats5.variables["x"][:]
ql = np.sum(stats5.variables["ql"][6,:,:,:],axis=0)
ql2 = np.sum(stats5.variables["ql"][24,:,:,:],axis=0)

X5,Y5 = np.meshgrid(x,y)


plt.figure(1)
plt.contourf(X5,Y5, ql, levels = 16, cmap= 'gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 3')
plt.suptitle('LWP at 100m')
plt.colorbar(cmap='gnuplot2', spacing = 'uniform')
plt.savefig('LWP1.png')

plt.figure(2)
plt.contourf(X5,Y5, ql2, levels= 16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 12')
plt.suptitle('LWP at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('LWP2.png')

#Rain Water Path
t = stats6.variables["time"][:]
z = stats6.variables["z"][:]
y = stats6.variables["y"][:]
x = stats6.variables["x"][:]
qr = np.sum(stats6.variables["qr"][6,:,:,:],axis=0)
qr2 = np.sum(stats6.variables["qr"][24,:,:,:],axis=0)

X6,Y6 = np.meshgrid(x,y)


plt.figure(9)
plt.contourf(X6,Y6, qr, levels = 16, cmap= 'gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 3')
plt.suptitle('RWP at 100m')
plt.colorbar(cmap='gnuplot2', spacing = 'uniform')
plt.savefig('RWP1.png')

plt.figure(10)
plt.contourf(X6,Y6, qr2, levels= 16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 12')
plt.suptitle('RWP at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('RWP2.png')

#Temperature Plot

t = stats4.variables["time"][:]
z = stats4.variables["z"][:]
y = stats4.variables["y"][:]
x = stats4.variables["x"][:]
T = stats4.variables["T"][6,4,:,:]
T2 = stats4.variables["T"][24,4,:,:]

X4,Y4 = np.meshgrid(x,y)

plt.figure(11)
plt.contourf(X4,Y4, T, levels = 16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 3')
plt.suptitle('Temperature at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('Temperature1.png')

plt.figure(12)
plt.contourf(X4,Y4,T2, levels=16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 12')
plt.suptitle('Temperature at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('Temperature2.png')

#Potential Temperature (Liquid Water)

t = stats7.variables["time"][:]
z = stats7.variables["z"][:]
y = stats7.variables["y"][:]
x = stats7.variables["x"][:]
thl = stats7.variables["thl"][6,4,:,:]
thl2 = stats7.variables["thl"][24,4,:,:]

X7,Y7 = np.meshgrid(x,y)

plt.figure(13)
plt.contourf(X7,Y7,thl, levels=16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 3')
plt.suptitle('Liquid Water Potential Temp at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('LWPotTemp1.png')

plt.figure(14)
plt.contourf(X7,Y7,thl2, levels=16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 12')
plt.suptitle('Liquid Water Potential Temp at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('LWPotTemp2.png')



#Wind Velocity Plots
#U-Velocity (X-Direction)
t = stats2.variables["time"][:]
z = stats2.variables["z"][:]
y = stats2.variables["y"][:]
xh = stats2.variables["xh"][:]
u = stats2.variables["u"][6,4,:,:]
u2 = stats2.variables["u"][24,4,:,:]

X2,Y2 = np.meshgrid(xh,y)


plt.figure(3)
plt.contourf(X2,Y2,u, levels = 16,cmap = 'gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 3')
plt.suptitle('U-Velocity at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('u-velocity1.png')

plt.figure(4)
plt.contourf(X2,Y2,u2, levels = 16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 12')
plt.suptitle('U-Velocity at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('u-velocity2.png')

#V-Velocity (Y-Direction)
t = stats3.variables["time"][:]
z = stats3.variables["z"][:]
yh = stats3.variables["yh"][:]
x = stats3.variables["x"][:]
v = stats3.variables["v"][6,4,:,:]
v2 = stats3.variables["v"][24,4,:,:]

X3,Y3 = np.meshgrid(x,yh)

plt.figure(5)
plt.contourf(X3,Y3, v , levels = 16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 3')
plt.suptitle('V-Velocity at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('v-velocity1.png')

plt.figure(6)
plt.contourf(X3,Y3, v2, levels = 16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 12')
plt.suptitle('V-Velocity at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('v-velocity2.png')

#W-Velocity (Z-Direction)
t = stats1.variables["time"][:]
zh = stats1.variables["zh"][:]
y = stats1.variables["y"][:]
x = stats1.variables["x"][:]
w = stats1.variables["w"][6,4,:,:]
w2 = stats1.variables["w"][24,4,:,:]

X1,Y1 = np.meshgrid(x,y)

plt.figure(7)
plt.contourf(X1,Y1,w, levels = 16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 3')
plt.suptitle('W-Velocity at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('w-velocity1.png')

plt.figure(8)
plt.contourf(X1,Y1,w2, levels=16, cmap='gnuplot2')
plt.xlabel('x-direction (m)')
plt.ylabel('y-direction (m)')
plt.title('Time = Hour 12')
plt.suptitle('W-Velocity at 100m')
plt.colorbar(cmap='gnuplot2',spacing='uniform')
plt.savefig('w-velocity2.png')

plt.show()


