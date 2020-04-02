"""
Garrett Witzke
Constrain CAO statistics 
"""

import numpy as np
import struct
import netCDF4
import math
from matplotlib import pylab as plt

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
nr_path  = default.groups["default"].variables["nr_path"][:]        #Number Density Rain Path 
rainflx = default.groups["thermo"].variables["rr"][:]           #Mean Surface Rain Rate
thlflux = default.groups["default"].variables["thl_w"][:,:]     #turbulent flux of theta_l
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
tket   = 0.5*(u2t + v2t + 0.5*(w2t[:,0:-1] + w2t[:,1::]))
    

#--------------------------------Caluclations-----------------------------------------
start = 0 
end = t.size

#Cloud base and height

bottomArr    = np.zeros(t.size)                        #create empty array to save cloud base values into
topArr       = np.zeros(t.size)                        #create empty array to save cloud top values into
mean_topArr  = np.zeros(t.size)                        #create empty arry to save mean cloud top values into
for tIdx in range(t.size):                             #for each timestep
    isCloud = (np.where(qlt[tIdx,:] > 10e-5))[0]       #Search for ql values bases on threshold
    if len(isCloud) == 0:                              #There are no ql values at time t=0
        bottomArr[tIdx] = 0                            #so set first number in heightArray = 0
        continue                                       #continue to next interation of for loop
    cb = isCloud[0]                                    #Cloud base values
    ct = isCloud[-1]                                   #Cloud top values  
    mean_ct = np.mean(isCloud)                  
    #print('base =',z[cb])
    #print('top =',z[ct])
    #print('mean =',mean_ct)
    bottomArr[tIdx] = z[cb]                            #Dumping values into the empty arrays
    topArr[tIdx]    = z[ct]                            #Dumping values into the empty arrays
    mean_topArr[tIdx]  = z[int(mean_ct)]                 #Dumping values into the empty arrays

#Rain Water Path

#Precipitation flux at the surface

#Cloud Cover


#SHF
SHF = cp * rhoh * thlflux

SHF_avg = np.zeros(t.size)

for tIdx in range(t.size):
    SHF_avg[tIdx] = np.mean(SHF[tIdx,:])


#LHF
LHF = rhoh * qtflux

LHF_avg = np.zeros(t.size)

for tIdx in range(t.size):
    LHF_avg[tIdx] = np.mean(LHF[tIdx,:])

#Entrainment velocity
#First smooth out the BLD line 
#then take derivative of BLD and subtract off mean w-velocitys at each corresponding height


#Plot a best-fit line to the BLD data
poly_degree = 2                         

coeffs = np.polyfit(t,zi, poly_degree)
poly_eqn = np.poly1d(coeffs)
y_hat = poly_eqn(t)

deriv = np.gradient(y_hat)


#Turbulent Kinetic Energy 

TKE = np.zeros(t.size)

for tIdx in range(t.size):
    TKE[tIdx] = np.mean(tket[tIdx,:])

#ql + qr

ql_rt = qlt +qrt



#----------------------------------Plotting Cloud Base and Height---------------------------
#                                           Figure 3
#                                        De Roode et Al.

f = 1
plt.figure(f)
plt.plot(t/3600,bottomArr, '--',label ='Cloud Base')
plt.plot(t/3600,topArr, label = 'Cloud Top')
plt.plot(t/3600,mean_topArr, label = 'Mean Cloud Top')
plt.title('Cloud Height vs Time')
plt.xlabel('Time (hrs)')
plt.ylabel('Height (m)')
plt.legend(loc='upper left')


#---------------------------------Plotting Time Series-------------------------------------
#                                       Figure 4
#                                   De Roode et Al.

#Cloud Cover
f += 1
plt.figure(f)
plt.plot(t/3600, ql_cover)
plt.title('Cloud Cover')
plt.xlabel('Time (hrs)')
plt.ylabel('Cloud Cover')


f += 1
plt.figure(f)
plt.plot(t/3600, LWPt, label = 'default.nc')
plt.title('Liquid Water Path')
plt.xlabel('Time (hrs)')
plt.ylabel('LWP')
plt.xlim(0,14.5)


f += 1
plt.figure(f)
plt.plot(t/3600,RWPt, label = 'qr_path')
plt.title('Rain Water Path')
plt.xlabel('Time (hrs)')





"""
f += 1
plt.figure(f)
plt.plot(t/3600,SHF_avg)
plt.title('Sensible Heat Flux')
plt.xlabel('Time (hrs)')
plt.ylabel('SHF')

f += 1
plt.figure(f)
plt.plot(t/3600, LHF_avg)
plt.title('Latent Heat Flux')
plt.xlabel('Time (hrs)')
plt.ylabel('LHF')

"""
f += 1
plt.figure(f)
plt.plot(t/3600,rainflx)
plt.title("Mean Surface Rain Rate")     #Precipitatio flux at the surface?
plt.xlabel("Time (hrs)")


f += 1
plt.figure(f)
plt.plot(t/3600, TKE)
plt.title('Turbulent Kinetic Energy')
plt.xlabel('Time (hrs)')
plt.ylabel('TKE')

"""
f += 1
plt.figure(f)
plt.plot(t/3600, W_ent)
plt.title('Entrainment Velocity')
plt.xlabel('Time (hrs)')
plt.ylabel('Entrainment Velocity')
"""

f += 1
plt.figure(f)
plt.plot(t/3600, zi, label = 'BLD')
plt.plot(t/3600, y_hat, label = "Regression Line")
plt.title('Boundary layer Depth')
plt.xlabel('Time (hrs)')
plt.ylabel('BL Depth (m)')
plt.legend(loc="upper left")



#----------------------Plotting Horizontal Mean profiles at t = 12hrs---------------
#                           t[144] --> hour 12 (43200 seconds)
#                                       Figure 5 
#                                   De Roode et Al.


fig, axs = plt.subplots(1,4, figsize = (12,6), sharey=True, constrained_layout=True)
fig.suptitle('Horizontal Mean Profiles of Instantaneous fields at t = 12 hr')
axs[0].plot(thl[144,:32],z[:32]/1000)
axs[0].set_xlabel(r'$\theta_l\/(K)$')
axs[0].set_ylabel(r'$Height\/(m)$')
axs[0].axis([270,290,0,3])
axs[1].plot(qtt[144,:32],z[:32]/1000)
axs[1].set_xlabel(r'$q_v + q_l + q_r\/(g/kg)$')
axs[1].axis([-1,4,0,3])
axs[2].plot(ql_rt[144,:32], z[:32]/1000)
axs[2].set_xlabel(r'$q_l + q_r\/(g/kg)$')
axs[2].axis([0,0.5,0,3])
axs[3].plot(ql_frac[144,:32], z[:32]/1000)
axs[3].set_xlabel(r'$Cloud\/Fraction$')
axs[3].axis([0,1,0,3])



#------------------------Comparing Horizontal Means at t=10minutes and t=10 hours------------
#                                               Figure 6
#                                           De Roode et Al.



#-------------------------------Contour Plots---------------------------------------------
#                               Figures 8 & 9
#                              De Roode et Al.

stats1  = netCDF4.Dataset('w.nc', 'r')
stats2  = netCDF4.Dataset('u.nc', 'r')
stats3  = netCDF4.Dataset('v.nc', 'r')
stats4  = netCDF4.Dataset('T.nc', 'r')
stats5  = netCDF4.Dataset('ql.nc', 'r')


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

CS1 = ax1.contourf(X2,Y2,u,levels=7,cmap='gist_rainbow',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
cbar1 = fig.colorbar(CS1, ax = ax1)
CS2 = ax2.contourf(X2,Y2,u2,levels=7, cmap='gist_rainbow')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
cbar2 = fig.colorbar(CS2, ax = ax2)


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

CS1 = ax1.contourf(X3,Y3,v,levels=7,cmap='gist_rainbow',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
cbar1 = fig.colorbar(CS1, ax = ax1)
CS2 = ax2.contourf(X3,Y3,v2,levels=7, cmap='gist_rainbow')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
cbar2 = fig.colorbar(CS2, ax = ax2)



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

CS1 = ax1.contourf(X1,Y1,w,levels=7,cmap='gist_rainbow',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
cbar1 = fig.colorbar(CS1, ax = ax1)
CS2 = ax2.contourf(X1,Y1,w2,levels=7, cmap='gist_rainbow')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
cbar2 = fig.colorbar(CS2, ax = ax2)

#Contour Plot pr T

t = stats4.variables["time"][:]
z = stats4.variables["z"][:]
y = stats4.variables["y"][:]
x = stats4.variables["x"][:]
T = stats4.variables["T"][6,5,:,:]
T2 = stats4.variables["T"][24,5,:,:]

X4,Y4 = np.meshgrid(x,y)

fig, (ax1,ax2) = plt.subplots(ncols=2, squeeze=True, sharey=True, figsize=(8,4))
fig.suptitle("T (K) @ 100m")

CS1 = ax1.contourf(X4,Y4,T,levels=7,cmap='gist_rainbow',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
cbar1 = fig.colorbar(CS1, ax = ax1)
CS2 = ax2.contourf(X4,Y4,T2,levels=7, cmap='gist_rainbow')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
cbar2 = fig.colorbar(CS2, ax = ax2)

#Contour Plot for LWP(ql_path) but for now its just ql
t = stats5.variables["time"][:]
z = stats5.variables["z"][:]
y = stats5.variables["y"][:]
x = stats5.variables["x"][:]
ql = stats5.variables["ql"][6,5,:,:]
ql2 = stats5.variables["ql"][24,5,:,:]

X5,Y5 = np.meshgrid(xh,y)

fig, (ax1,ax2) = plt.subplots(ncols=2, squeeze=True, sharey=True, figsize=(8,4))
fig.suptitle("ql @ 100m")

CS1 = ax1.contourf(X5,Y5,ql,levels=7,cmap='gist_rainbow',)
ax1.set_title('Time = Hour 3')
ax1.set_xlabel('x-direction (m)')
cbar1 = fig.colorbar(CS1, ax = ax1)
CS2 = ax2.contourf(X5,Y5,ql2,levels=7, cmap='gist_rainbow')
ax2.set_title('Time = Hour 12')
ax2.set_ylabel('y-direction (m)')
cbar2 = fig.colorbar(CS2, ax = ax2)







plt.show()





  