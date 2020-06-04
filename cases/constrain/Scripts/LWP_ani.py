"""
Garrett Witzke
Constrain CAO statistics 
Contour Plots
"""

import numpy as np
import netCDF4
from matplotlib import pylab as plt
from matplotlib import animation as animation
from matplotlib.animation import PillowWriter


#----------------------------Constants----------------------------------------------
Rd = 287.0                #Universal Gas constant for dry air (J/K mol)
cp = 1006.0               #Specific Heat of Dry Air (J/kg K)
Lv = 2.26e6               #Latent Heat of Vaporization (J/kg)
p0 = 100900               #Surface pressure (Pa)
f  = 1.3e-4               #1/seconds Coriolis Parameter at a latitude of 63
Nc = 50                   #1/cm^3 Cloud Droplet conc


#----------------------------import results-----------------------------------------
stats5  = netCDF4.Dataset('ql.nc', 'r')


#-------------------------------Contour Plots---------------------------------------------
#                               Figures 8 & 9
#                              De Roode et Al.
#                        t[6]= hour 3, t[24]= hour 12

#Liquid Water Path
fig = plt.figure()

t = stats5.variables["time"][:]
z = stats5.variables["z"][:]
y = stats5.variables["y"][:]
x = stats5.variables["x"][:]
ql_field = stats5.variables["ql"][:,:,:,:]

X5,Y5 = np.meshgrid(x,y)
dz = z[-1]/z.size
ql_path = np.zeros((x.size,y.size))

ims=[]
for t in range(t.size):
    ql_path = np.sum(ql_field[t,:,:,:] * dz ,axis=0)
    ql_path[ql_path==0]=np.nan
    im = plt.imshow(ql_path, cmap='gnuplot2', animated=True, origin='lower', extent=(x[0],x[-1],y[0],y[-1]))
    plt.title(r'$Liquid\/Water\/Path\/(kg/m^2)$')
    plt.xlabel('x-direction (m)')
    plt.ylabel('y-direction (m)')


    ims.append([im])

ani = animation.ArtistAnimation(fig,ims, interval = 50, blit = True, repeat_delay = 500)

writer = PillowWriter(fps=1)
ani.save('LWP.gif', writer=writer)

plt.show()