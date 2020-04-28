# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 11:13:44 2018

@author: krueg

Using years with anomalously high surface mass balance (as determined
previously), determine sea surface temperature anomalies as compared to the
long-term average (1850-2010).
"""

# Import modules
import time
start_time = time.time()

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap


# Load file
def load_file():
    filename = 'sst.mon.mean.nc'
    return Dataset(filename,"r")

nc = load_file()


# Initialize variables
def temp_var():
    return load_file().variables['sst']
    # Shape #(2016, 180, 360)[months since 1850-01-01, lat, long]
    
temp = temp_var()


# Earth-sized array initialization
newshape = [temp.shape[1],temp.shape[2]] # (180, 360)


# Detrend SST values
def detrend_sst():
    return signal.detrend(temp, axis=0, type='linear') 
    # Shape (2016, 180, 360)

d_sst = detrend_sst()

# Calculate long-term mean of SST (1850-2010)
def mean_sst():
    """Calculate long-term mean of SST (1850-2010)"""
    
    avg = np.zeros(newshape)   # (180, 360)
    
    for i in range(0, 1932):  # 1850-2010 by month
    # Read temperatures into each grid cell for all months additive
        avg += d_sst[i, :, :]
    # Divide each grid value by length of observations (1932 months = 161 years)
    avg = avg / 1932
    # avg.shape #Shape is (180, 360)
    
    # Remove nan cells from dataset
    for x in range(0, 180):
        for y in range(0, 360):
            if avg[x, y] > 100.: # Missing values are 1e20
                avg[x, y] = np.nan
                
    return avg

avg = mean_sst()
                    

# Find mean SST for anomalously high SMB years 
def high_sst():
    """Calculate mean SSTs during high SMB years."""
    high_years = [1866, 
                  1900, 
                  1901, 
                  1902, 
                  1984, 
                  1985, 
                  1988, 
                  1992, 
                  1998, 
                  1999, 
                  2000]
    base_year = 1850
    years = len(high_years)
    high_years_upper = []
    high_years_lower = []
    for i in range(0, years):
        high_years_upper.append((high_years[i] - base_year) * 12)
        high_years_lower.append(high_years_upper[i] + 12)
    
    high = np.zeros(newshape)   # (180, 360)
    for i in range(0, years):
        for j in range(high_years_upper[i], high_years_lower[i]):
            high += d_sst[j, :, :]
    # Divide each value by length of observations (years * 12 months)
    high = high / (years * 12)
    
    # Remove nan cells 
    for x in range(0, 180):
        for y in range(0, 360):
            if high[x, y] > 100.: # Missing values are 1e20
                high[x, y] = np.nan
    
    return high

high = high_sst()
     
          
#Find SST deviations of anomalously high SMB years from long-term average SST
def sst_devs():         
    sst_dev = high[:, :] - avg[:, :]
    print(sst_dev)
    
    return sst_dev

sst_dev = sst_devs()


# Perform student's t-test (two-tailed 95%)
# 10 degrees of freedom (11 observation years)
def sig_devs():
    mean = np.nanmean(sst_dev)
    stddev = np.nanstd(sst_dev)
    n = 11
    mu_plus = mean + 2.228 * (stddev / (n**0.5))
    mu_minus = mean - 2.228 * (stddev / (n**0.5))
    print(mean, stddev, n, mu_plus, mu_minus)
    sig = np.zeros(newshape)
    for x in range(0, 180):
        for y in range(0, 360):
            if sst_dev[x, y] > mu_plus:
                sig[x, y] = 1.00000000001
            elif sst_dev[x, y] < mu_minus:
                sig[x, y] = -1.00000000001
            else: 
                sig[x, y] = 0 
                
    return sig

sig = sig_devs()


# Add core locations to plot
import csv

def core_lons():
    cores = open('group_wais_cores.csv', 'r')
    reader = csv.reader(cores)
    all_values = []
    for row in reader:
        all_values.append(row)
    cores.close()
    all_values = np.array(all_values)
    
    core_lons = []
    for i in range(1, 40):
        name = all_values[0, i]
        lon = float(all_values[2, i])
        core_lons.append(lon)
    
    return core_lons

core_lon = core_lons()

def core_lats():
    cores = open('group_wais_cores.csv', 'r')
    reader = csv.reader(cores)
    all_values = []
    for row in reader:
        all_values.append(row)
    cores.close()
    all_values = np.array(all_values)
    
    core_lats = []
    for i in range(1, 40):
        name = all_values[0, i]
        lat = float(all_values[1, i])
        core_lats.append(lat)
    
    return core_lats

core_lat = core_lats()


# Plot SST deviations of anomalously high SMB years from long-term average SST
from textwrap import wrap

def plot_results():
    lats = nc.variables['lat'][:]
    lons = nc.variables['lon'][:]
    print(lons.min(),lons.max())
    print(lats.min(),lats.max())
    
    variable = sst_dev[:,:] # Plotting anomalously high SMB year SST deviations
    print(np.nanmax(sst_dev), np.nanmin(sst_dev))
    # (0.38,-0.49 using detrended data)
    
    lonmin, lonmax, latmin, latmax = lons.min(), lons.max(), lats.min(), lats.max()
    
    m = Basemap(llcrnrlat=latmin, urcrnrlat=45., llcrnrlon=20., urcrnrlon=300., 
                resolution = 'h', rsphere=(6367470.00, 6367470.00),
                projection = 'cyl')
    
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)
    
    cs = m.pcolormesh(xi, yi, variable, cmap='bwr', vmin=-0.8, vmax=0.8)
    m.drawparallels(np.arange(-80., 90., 10.), labels=[1, 0, 0, 0], fontsize=16)
    m.drawmeridians(np.arange(-180., 180., 20.), labels=[0, 0, 0, 1], fontsize=16)
    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color='lightgrey')
    m.contour(xi, yi, sig[:, :], [-1., 1.], colors='k', linewidths=3)
    cbar = m.colorbar(cs, location = 'bottom', pad = '10%')
    c = "SST Anomalies Relative to 1850-2010 Mean (Degrees C)"
    cbar.set_label('\n'.join(wrap(c, 45)), size=24)
    title = "SST Anomalies in High SMB Years Relative to 1850-2010 SST Mean"  
    plt.title('\n'.join(wrap(title, 45)), size=32) 
    
    x, y = m(core_lon, core_lat)
    m.scatter(x, y, zorder=10, latlon=True, color='k')
    
    return plt.show()


plot_results()

print("My program took", time.time() - start_time, "seconds to run")