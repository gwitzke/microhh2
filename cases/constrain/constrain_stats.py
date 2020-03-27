"""
Garrett Witzke
Constrain CAO statistics 
"""

import numpy as np
import struct
import netCDF4
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
stats_ql = netCDF4.Dataset("constrain_ql_0000000.nc", 'r')

#Only grabbing necessary data

#Assigning basic variable names
t    = default.variables["time"][:]                              #Time
z    = default.variables["z"][:]                                 #Height 
zh   = default.variables["zh"][:38]                              #Half Height
zi   = default.groups["thermo"].variables["zi"][:]               #Boundary layer depth
zi1  = stats_ql.groups["thermo"].variables["zi"][:]  
rho  = default.groups["thermo"].variables["rho"][:,:]            #Density
rhoh = default.groups["thermo"].variables["rhoh"][:,:]           #Half level Density
Tt   = default.groups["thermo"].variables["T"][:,:]              #Absolute Temperature

#Assigning fractional area contained in mask
areat  = default.groups["default"].variables["area"][:,:]
areaht = default.groups["default"].variables["areah"][:,:]

#Assigning Time Series Variable names
LWPt  = default.groups["thermo"].variables["ql_path"][:]*1000   #Liquid Water Path units of g/m^2
LWPt1 = stats_ql.groups["thermo"].variables["ql_path"][:]*1000
thlflux = default.groups["default"].variables["thl_w"][:,:]     #turbulent flux of theta_l
qtflux = default.groups["default"].variables["qt_w"][:,:]       #Turbulent flux of qt                      
    
#Assigning Thermodynamic variable names
thl  = default.groups["default"].variables["thl"][:,:]        #Theta_l
qtt  = default.groups["default"].variables["qt"][:,:]*1000    #qt units of g/kg
qlt  = default.groups["thermo"].variables["ql"][:,:]*1000     #ql units of g/kg
qlt1 = stats_ql.groups["thermo"].variables["ql"][:,:]*1000 
ql_cover = default.groups["thermo"].variables["ql_cover"][:]
ql_cover1 = stats_ql.groups["thermo"].variables["ql_cover"][:]
qvt  = qtt - qlt
ql_frac = default.groups["thermo"].variables["ql_frac"][:,:]
ql_frac1 = stats_ql.groups["thermo"].variables["ql_frac"][:,:]


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

bottomArr = np.zeros(t.size)                        #create empty array to save cloud base values into
topArr    = np.zeros(t.size)                        #create empty array to save cloud top values into
mean_top  = np.zeros(t.size)                        #create empty arry to save mean cloud top values into
for tIdx in range(t.size):                          #for each timestep
    isCloud = (np.where(qlt[tIdx,:] > 10e-5))[0]    #Search for ql values bases on threshold
    if len(isCloud) == 0:                           #There are no ql values at time t=0
        bottomArr[tIdx] = 0                         #so set first number in heightArray = 0
        continue                                    #continue to next interation of for loop
    cb = isCloud[0]                                 #Cloud base values
    ct = isCloud[-1]                                #Cloud top values
    mean_ct = np.mean(isCloud)                      #Mean Cloud top values
    #print(z[cb])
    #print'cloud top height =',z[ct]
    #print'Time steps for heights where clouds are located',isCloud
    #print(mean_top)
    bottomArr[tIdx] = z[cb]                         #
    topArr[tIdx]    = z[ct]                         # 
    mean_top[tIdx]  = z[mean_ct]                    #

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
W_h = np.zeros(t.size)
for Idx in range(t.size):
    W_h[Idx] = np.mean(wt[Idx,:])
#print(W_h)
#print(wt)
W_ent = np.gradient(zi) - W_h

#Turbulent Kinetic Energy 

TKE = np.zeros(t.size)

for tIdx in range(t.size):
    TKE[tIdx] = np.mean(tket[tIdx,:])



#----------------------------------Plotting Cloud Base and Height---------------------------
#                                           Figure 3
#                                        De Roode et Al.

f = 1
plt.figure(f)
plt.plot(t/3600,bottomArr, '--',label ='Cloud Base')
plt.plot(t/3600,topArr, label = 'Cloud Top')
plt.plot(t/3600,mean_top, label = 'Mean Cloud Top')
plt.title('Cloud Height vs Time')
plt.xlabel('Time (hrs)')
plt.ylabel('Height (m)')
plt.legend(loc='upper left')


#---------------------------------Plotting Time Series-------------------------------------
#                                       Figure 4
#                                   De Roode et Al.



f += 1
plt.figure(f)
plt.plot(t/3600, qlt1)
plt.xlabel('Time (hrs)')
plt.ylabel('ql')







#Cloud Cover
f += 1
plt.figure(f)
plt.plot(t/3600, ql_cover, label = 'default.nc')
plt.plot(t/3600,ql_cover1, label = 'ql.nc', ls = '--')
plt.title('Cloud Cover')
plt.xlabel('Time (hrs)')
plt.ylabel('Cloud Cover')
plt.legend(loc='lower right')


f += 1
plt.figure(f)
plt.plot(t/3600, LWPt, label = 'default.nc')
plt.plot(t/3600, LWPt1, label = 'ql.nc', ls = '--')
plt.title('Liquid Water Path')
plt.xlabel('Time (hrs)')
plt.ylabel('LWP')
plt.xlim(0,14.5)
plt.legend(loc='lower right')

#Rain Water Path---After warm micro added

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

#Precipitation flux at surface

f += 1
plt.figure(f)
plt.plot(t/3600, TKE)
plt.title('Turbulent Kinetic Energy')
plt.xlabel('Time (hrs)')
plt.ylabel('TKE')


f += 1
plt.figure(f)
plt.plot(t/3600, W_ent)
plt.title('Entrainment Velocity')
plt.xlabel('Time (hrs)')
plt.ylabel('Entrainment Velocity')
"""

f += 1
plt.figure(f)
plt.plot(t/3600, zi, label = 'default.nc')
plt.plot(t/3600, zi1, label = 'ql.nc', ls = '--')
plt.title('Boundary layer Depth')
plt.xlabel('Time (hrs)')
plt.ylabel('BL Depth (m)')
plt.legend(loc='lower right')


#----------------------Plotting Horizontal Mean profiles at t = 12hrs---------------
#                           t[144] --> hour 12 (43200 seconds)
#                                       Figure 5 
#                                   De Roode et Al.


fig, axs = plt.subplots(1,4, figsize = (12,6), sharey=True, constrained_layout=True)
fig.suptitle('Horizontal Mean Profiles of Instantaneous fields at t = 12 hr')
axs[0].plot(thl[144],z[:]/1000)
axs[0].set_xlabel(r'$\theta_l\/(K)$')
axs[0].set_ylabel(r'$Height\/(m)$')
axs[0].axis([270,290,0,5])
axs[1].plot(qtt[144,:],z[:]/1000)
axs[1].set_xlabel(r'$q_v + q_l + q_r\/(g/kg)$')
axs[1].axis([0,4,0,5])
axs[2].plot()
axs[2].set_xlabel(r'$q_l + q_r\/(g/kg)$')
axs[3].plot()
axs[3].set_xlabel(r'$Cloud\/Fraction$')



#------------------------Comparing Horizontal Means at t=10minutes and t=10 hours------------
#                                               Figure 6
#                                           De Roode et Al.



#-------------------------------Contour Plots---------------------------------------------
#                               Figures 8 & 9
#                              De Roode et Al.






plt.show()





  