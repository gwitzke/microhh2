
import numpy as np
import struct
import netCDF4
from matplotlib import pylab as plt


default = netCDF4.Dataset("constrain_default_0000000.nc", 'r') 
qlt  = default.groups["thermo"].variables["ql"][:,:]*1000     #ql units of g/kg
t    = default.variables["time"][:]   
z    = default.variables["z"][:]  


#Cloud base and height

bottomArr = np.zeros(t.size)                        #create empty array to save cloud base values into
topArr    = np.zeros(t.size)                        #create empty array to save cloud top values into
mean_top  = np.zeros(t.size)                        #Create empty array to save mean cloud top values into
for tIdx in range(t.size):                          #for each timestep
    isCloud = (np.where(qlt[tIdx,:] > 10e-5))[0]    #Search for ql values bases on threshold
    if len(isCloud) == 0:                           #There are no ql values at time t=0
        bottomArr[tIdx] = 0                         #so set first number in heightArray = 0
        continue                                    #continue to next interation of for loop
    cb = isCloud[0]                                 #Cloud base values
    ct = isCloud[-1]                                #Cloud top values
    mean_ct = np.mean(isCloud)
    #print(z[cb])
    print'cloud top height =',z[ct]
    print'Time steps for heights where clouds are located',isCloud
    print(mean_top)
    bottomArr[tIdx] = z[cb]                         #
    topArr[tIdx]    = z[ct]                         # 
    mean_top[tIdx]  = z[mean_ct]


f = 1
plt.figure(f)
plt.plot(t/3600,bottomArr, '--',label ='Cloud Base')
plt.plot(t/3600,topArr, label = 'Cloud Top')
plt.plot(t/3600,mean_top, label='Mean Cloud Top')
plt.title('Cloud Height vs Time')
plt.xlabel('Time (hrs)')
plt.ylabel('Height (m)')
plt.legend(loc='upper left')
plt.show()