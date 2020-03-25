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
stats = netCDF4.Dataset("constrain_default_0000000.nc", 'r') 

#Only grabbing necessary data

#Assigning basic variable names
t    = stats.variables["time"][:]                              #Time
z    = stats.variables["z"][:]                                 #Height 
zh   = stats.variables["zh"][:38]                              #Half Height
zi   = stats.groups["thermo"].variables["zi"][:]               #Boundary layer depth
rho  = stats.groups["thermo"].variables["rho"][:,:]            #Density
rhoh = stats.groups["thermo"].variables["rhoh"][:,:]           #Half level Density
Tt   = stats.groups["thermo"].variables["T"][:,:]              #Absolute Temperature

#Assigning fractional area contained in mask
areat  = stats.groups["default"].variables["area"][:,:]
areaht = stats.groups["default"].variables["areah"][:,:]

#Assigning Time Series Variable names
LWPt  = stats.groups["thermo"].variables["ql_path"][:]*1000   #Liquid Water Path units of g/m^2
thlflux = stats.groups["default"].variables["thl_w"][:,:]     #turbulent flux of theta_l
qtflux = stats.groups["default"].variables["qt_w"][:,:]       #Turbulent flux of qt                      
    
#Assigning Thermodynamic variable names
thl  = stats.groups["default"].variables["thl"][:,:]        #Theta_l
qtt  = stats.groups["default"].variables["qt"][:,:]*1000    #qt units of g/kg
qlt  = stats.groups["thermo"].variables["ql"][:,:]*1000     #ql units of g/kg
qvt  = qtt - qlt
ql_frac = stats.groups["thermo"].variables["ql_frac"][:,:]


#Assigning Wind Component variable names
wt     = stats.groups["default"].variables["w"][:,:]        #Vertical Velocity (z-direction)
ut     = stats.groups["default"].variables["u"][:,:]        #U direction velocity (x- component)
u2t    = stats.groups["default"].variables["u_2"][:,:]      
vt     = stats.groups["default"].variables["v"][:,:]        #V direction velocity (y -component)
v2t    = stats.groups["default"].variables["v_2"][:,:]
w2t    = stats.groups["default"].variables["w_2"][:,:]
ufluxt = stats.groups["default"].variables["u_w"][:,:]      #Turblulent flux of u velocity
vfluxt = stats.groups["default"].variables["v_w"][:,:]      #Turbulent flux of v velocity
Ufluxt = ufluxt + vfluxt
tket   = 0.5*(u2t + v2t + 0.5*(w2t[:,0:-1] + w2t[:,1::]))
    

#--------------------------------Caluclations-----------------------------------------
start = 0 
end = t.size

#Cloud base and height

bottomArr = np.zeros(t.size)                        #create empty array to save cloud base values into
topArr    = np.zeros(t.size)                        #create empty array to save cloud top values into
for tIdx in range(t.size):                          #for each timestep
    isCloud = (np.where(qlt[tIdx,:] > 10e-5))[0]    #Search for ql values bases on threshold
    if len(isCloud) == 0:                           #There are no ql values at time t=0
        bottomArr[tIdx] = 0                         #so set first number in heightArray = 0
        continue                                    #continue to next interation of for loop
    cb = isCloud[0]                                 #Cloud base values
    ct = isCloud[-1]                                #Cloud top values
    #print(z[cb])
    #print(z[ct])
    bottomArr[tIdx] = z[cb]                         #
    topArr[tIdx]    = z[ct]                         # 

#Rain Water Path

#Precipitation flux at the surface

#Cloud Fraction

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


#Cloud Coverage and fraction

ql_frac_avg = np.zeros(t.size)
ql_avg = np.zeros(t.size)

for tIdx in range(t.size):
    ql_frac_avg[tIdx] = np.mean(ql_frac[tIdx,:])
    ql_avg[tIdx] = np.mean(qlt[tIdx,:])




#----------------------------------Plotting Cloud Base and Height---------------------------
#                                           Figure 3
#                                        De Roode et Al.

f = 1
plt.figure(f)
plt.plot(t/3600,bottomArr, '--',label ='Cloud Base')
plt.plot(t/3600,topArr, label = 'Cloud Top')
plt.title('Cloud Height vs Time')
plt.xlabel('Time (hrs)')
plt.ylabel('Height (m)')
plt.legend(loc='upper left')


#---------------------------------Plotting Time Series-------------------------------------
#                                       Figure 4
#                                   De Roode et Al.

f += 1
plt.figure(f)
plt.plot(t/3600, ql_frac_avg)
plt.title('Cloud Fraction')


f += 1
plt.figure(f)
plt.plot(t/3600, ql_avg)
plt.title('Cloud Covere')




f += 1
plt.figure(f)
plt.plot(t/3600, LWPt)
plt.title('Liquid Water Path')
plt.xlabel('Time (hrs)')
plt.ylabel('LWP')
plt.xlim(0,14.5)

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





  