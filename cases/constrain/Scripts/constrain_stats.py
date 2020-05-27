"""
Garrett Witzke
Constrain CAO statistics 
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
g  = 9.81                 #Gravitational Acceleration


#----------------------------import results-----------------------------------------
default = netCDF4.Dataset("constrain.default.0000000.nc", 'r') 



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
LWPt     = default.groups["thermo"].variables["ql_path"][:]*1000    #Liquid Water Path units of g/m^2
RWPt     = default.groups["thermo"].variables["qr_path"][:]*1000   #Rain Water Specific Humidity Path g/m^2
nr_path  = default.groups["thermo"].variables["nr_path"][:]        #Number Density Rain Path 
rainflux = default.groups["thermo"].variables["rr"][:]              #Mean Surface Rain Rate
thlflux  = default.groups["thermo"].variables["thl_w"][:,:]        #turbulent flux of theta_l
qtflux   = default.groups["thermo"].variables["qt_w"][:,:]         #Turbulent flux of qt                      
    
#Assigning Thermodynamic variable names
thl  = default.groups["thermo"].variables["thl"][:,:]        #Theta_l
qtt  = default.groups["thermo"].variables["qt"][:,:]*1000    #qt units of g/kg
qlt  = default.groups["thermo"].variables["ql"][:,:]*1000     #ql units of g/kg
qrt  = default.groups["thermo"].variables["qr"][:,:]*1000    #qr units of g/kg
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

#Assigning Radiation Component variables names
sflux  = default.groups["radiation"].variables["sflx"][:,:]    #Total Shortwave Radiative Flux
lflux  = default.groups["radiation"].variables["lflx"][:,:]    #Total Longwave Radiative Flux
    

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
    mean_topArr[tIdx]  = z[int(mean_ct)]               #Dumping values into the empty arrays


#Precipitation flux at the surface
    #Need to convert units to W/m^2
Rainy_CloudArr = np.zeros(t.size)
no_rainflx  = np.zeros(t.size)

for tIdx in range(t.size):
    RainyCloud = (np.where(qrt[tIdx,:] > 1e-4))[0]
    if len(RainyCloud) == 0:
        no_rainflx[tIdx] = 0
        continue
    rct = RainyCloud[-1]
    mean_rct = np.mean(RainyCloud)
    Rainy_CloudArr[tIdx] = z[int(mean_rct)]

sfc_rainflx = rainflux * g * Rainy_CloudArr


#Cloud Cover
    #ql_cover

#SHF

SHF = cp * rhoh * thlflux
SHF_avg = np.zeros(t.size)
SHF_max = np.zeros(t.size)
for tIdx in range(t.size):
    SHF_avg[tIdx] = - np.mean(SHF[tIdx,:])
    SHF_max[tIdx] = np.max(SHF[tIdx,:])

#LHF
    
LHF = Lv * rhoh * qtflux
LHF_avg = np.zeros(t.size)
LHF_max = np.zeros(t.size)
for tIdx in range(t.size):
    LHF_avg[tIdx] = np.mean(LHF[tIdx,:])
    LHF_max[tIdx] = np.max(LHF[tIdx,:])

#Entrainment velocity
#First smooth out the BLD line 
#then take derivative of BLD and subtract off mean w-velocitys at each corresponding height


        #Plot a best-fit line to the BLD data
poly_degree = 2                         

coeffs = np.polyfit(t,zi, poly_degree)
poly_eqn = np.poly1d(coeffs)
y_hat = poly_eqn(t)

deriv = np.gradient(y_hat,t)
W_ent = np.zeros(t.size)

for tIdx in range(t.size):
    for zIdx in range(zh.size):
        W_ent[tIdx] = deriv[tIdx] - wt[tIdx,zIdx]

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
plt.axis([0,14.5,0,3000])


#---------------------------------Plotting Time Series-------------------------------------
#                                       Figure 4
#                                   De Roode et Al.

#Subplotting them all together
fig, axs = plt.subplots(7,1, figsize=(6,10), sharex=True, tight_layout=True)
fig.suptitle("Figure 4 De Roode")
axs[0].plot(t/3600, ql_cover)
axs[0].set_ylabel('Cloud Cover')
axs[1].plot(t/3600, LWPt)
axs[1].set_ylabel(r'$LWP\/(g/m^2)$')
axs[2].plot(t/3600, RWPt)
axs[2].set_ylabel(r'$RWP\/(g/m^2)$')
axs[3].plot(t/3600, SHF_max)
axs[3].set_ylabel(r'$SHF\/(W/m^2)$')
axs[4].plot(t/3600, LHF_max)
axs[4].set_ylabel(r'$LHF\/(W/m^2)$')
axs[5].plot(t/3600, sfc_rainflx)
axs[5].set_ylabel(r'$sfc\/prec\/(W/m^2)$')
axs[6].plot(t/3600, deriv)
axs[6].set_ylabel(r'$W_e\/(m/s)$')
axs[6].set_xlabel("Time (hrs)")


"""
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
plt.ylabel(r'$LWP\/(g/m^2)$')
plt.xlim(0,14.5)


f += 1
plt.figure(f)
plt.plot(t/3600,RWPt, label = 'qr_path')
plt.title('Rain Water Path')
plt.xlabel('Time (hrs)')
plt.ylabel(r'$RWP\/(g/m^2)$')

f += 1
plt.figure(f)
plt.plot(t/3600,SHF_avg, label = 'Average')
plt.plot(t/3600,SHF_max, label = 'Max')
plt.title('Sensible Heat Flux')
plt.xlabel('Time (hrs)')
plt.ylabel(r'$SHF\/(W/m^2)$')
plt.legend(loc='upper right')

f += 1
plt.figure(f)
plt.plot(t/3600, LHF_avg, label = 'Average')
plt.plot(t/3600, LHF_max, label = 'Max')
plt.title('Latent Heat Flux')
plt.xlabel('Time (hrs)')
plt.ylabel(r'$LHF\/(W/m^2)$')
plt.legend(loc='upper right')


f += 1
plt.figure(f)
plt.plot(t/3600, sfc_rainflx, label = 'sfcflux')
#Splt.plot(t/3600, rainflux,     label = 'rr')
plt.title("Precipitation Flux at the Surface")  
plt.suptitle('Mean Surface Rain Rate')  
plt.xlabel("Time (hrs)")
plt.ylabel(r'$(W/m^2)$')
plt.legend(loc='upper left')



#f += 1
#plt.figure(f)
#plt.plot(t/3600, TKE)
#plt.title('Turbulent Kinetic Energy')
#plt.xlabel('Time (hrs)')
#plt.ylabel('TKE')


f += 1
plt.figure(f)
plt.plot(t/3600, deriv)
plt.title('Entrainment Velocity')
plt.xlabel('Time (hrs)')
plt.ylabel(r'$W_e\/(m/s)$')
plt.yticks((0.00,0.05,0.10))
plt.axis([0,14.5,0,0.10])

"""

plt.figure()
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
axs[0].plot(thl[23,:],z[:]/1000)
axs[0].set_xlabel(r'$\theta_l\/(K)$')
axs[0].set_ylabel(r'$Height\/(km)$')
axs[0].axis([270,290,0,3])
axs[1].plot(qtt[23,:],z[:]/1000)
axs[1].set_xlabel(r'$q_v + q_l + q_r\/(g/kg)$')
axs[1].axis([-1,4,0,3])
axs[2].plot(ql_rt[23,:], z[:]/1000)
axs[2].set_xlabel(r'$q_l + q_r\/(g/kg)$')
axs[2].axis([0,0.5,0,3])
axs[3].plot(ql_frac[23,:], z[:]/1000)
axs[3].set_xlabel(r'$Cloud\/Fraction$')
axs[3].axis([0,1,0,3])


#---------------------Plotting Radiation Terms from [Radiation]---------------------











plt.show()





  
