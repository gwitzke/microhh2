"""
Garrett Witzke
Constrain CAO statistics 
"""

import numpy as np
import struct
import netCDF4
from pylab import *

nt = 12
plotens = False

def plotstats(name,line):
    stats = netCDF4.Dataset("constrain_{}_0000000.nc".format(name), 'r') #import results

    #Assigning basic variable names
    t   = stats.variables["time"][:]
    z   = stats.variables["z"][:]
    zh  = stats.variables["zh"][:]

    #Assigning fractional area contained in mask
    areat  = stats.groups["default"].variables["area"][:,:]
    areaht = stats.groups["default"].variables["areah"][:,:]

    #Assigning Time Series Variable names
    LWPt = stats.groups["thermo"].variables["ql_path"][:]
    #RWP(rain water path) =
    SHFt = stats.groups["default"].variables["thl_diff"][:,:]
    LHFt = stats.groups["default"].variables["qt_diff"][:,:]
    #SPF(surface precipatation flux)
    #Entrainment

    #Assigning Thermodynamic variable names
    thl  = stats.groups["default"].variables["thl"][:,:]
    qtt  = stats.groups["default"].variables["qt"][:,:]*1000
    qlt  = stats.groups["thermo"].variables["ql"][:,:]*1000

    
    #Assigning Wind Component variable names
    u  = stats.groups["default"].variables["u"][:,:]
    u2 = stats.groups["default"].variables["u_2"][:,:]
    v  = stats.groups["default"].variables["v"][:,:]
    v2 = stats.groups["default"].variables["v_2"][:,:]

    

    #Caluclations
    end   = t.size
    start = t.size - nt

    #Time Series Calculations
    Clouds = np.mean(areat[start:end,:], 0)
    #SHF    = np.sum(areat[start:end,:]*SHFt[start:end,:], 0) / np.sum(areat[start:end,:], 0)
    #LHF    = np.sum(areat[start:end,:]*LHFt[start:end,:], 0) / np.sum(areat[start:end,:], 0)
    LWP    = np.sum(areat[start:end]*LWPt[start:end], 0) / np.sum(areat[start:end], 0)



    #Thermodynamic Profile Calculations
    th_l = np.sum(areat[start:end,:]*thl[start:end,:], 0) / np.sum(areat[start:end,:], 0)
    qt   = np.sum(areat[start:end,:]*qtt[start:end,:], 0) / np.sum(areat[start:end,:], 0)
    ql   = np.sum(areat[start:end,:]*qlt[start:end,:], 0) / np.sum(areat[start:end,:], 0)

    #Plotting
    f = 1

    #Plots of Time Series
    figure(f)
    if(plotens):
        for n in range(start,end):
            plot(areat[n,:],z)
    plot(Clouds, z, line, label=name)
    title("Cloud Coverage")

    f += 1
    figure(f)
    if(plotens):
        for n in range(start,end):
            plot(LWPt[:],t)
    plot(LWP, t, line, label=name)
    title("Liqiud Water Path")


    #Plots of Thermodynamic profiles
    f += 1
    figure(f)
    if(plotens):
        for n in range(start, end):
            plot(thl[n,:], z)
    plot(th_l, z, line, label=name)
    title("Liquid Potential Temperature")

plotstats('default', 'k-')
show()










  