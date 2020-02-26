"""
Garrett Witzke
Constrain CAO statistics 
"""

import numpy as np
import struct
import netCDF4
import matplotlib as plt

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
    Cloud_cover
    LPW
    RWP
    SHF
    LHF
    SPF
    Entrain

    #Assigning Thermodynamic variable names
    th_l 
    qt 
    
    #Assigning Wind Component variable names
    u 
    v
    w

    #Caluclations
    end = t.size
    start = t.size -nt

    #Time Series Calculations


    #Thermodynamic Profile Calculations

    #Plotting
    f = 1

    #Plots of Time Series

    #Plots of Thermodynamic profiles

plotstats()
plt.show()










  