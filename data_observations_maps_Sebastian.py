#!/usr/bin/env python
# coding: utf-8

# In[ ]:

def plot_map(data, lat, lon, varname, cbar_label):
    lat_min = -60  # Minimum latitude for South America
    lat_max = 15   # Maximum latitude for South America
    lon_min = -90  # Minimum longitude for South America
    lon_max = -30  # Maximum longitude for South America

    data_avg=data.squeeze()

    m = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max, resolution='l') # Create the Basemap object for South America
    lon_mesh, lat_mesh = np.meshgrid(lon, lat)   # Create a meshgrid of latitudes and longitudes
    x, y = m(lon_mesh, lat_mesh)                 # Convert latitudes and longitudes to x, y coordinates
    plt.figure(figsize=(6, 6))
    
    if varname == 'TC': m.pcolormesh(x, y, data_avg, cmap='YlGn', latlon=True)
    #elif varname == 'PR': m.pcolormesh(x, y, np.where((data_avg>0)&(data_avg<=4000), data_avg, np.nan), cmap='YlGnBu', latlon=True)
    elif varname == 'PR': m.pcolormesh(x, y, data_avg, cmap='YlGnBu', vmax=4000, latlon=True)
    elif varname == 'TM': m.pcolormesh(x, y, data_avg, cmap='RdYlBu_r', vmin=15, latlon=True)
    elif varname == 'SW': m.pcolormesh(x, y, data_avg, cmap='YlOrBr', vmax=250, latlon=True)
    elif varname == 'LW': m.pcolormesh(x, y, data_avg, cmap='YlOrBr', vmin=300, latlon=True)
    elif varname == 'SD': m.pcolormesh(x, y, data_avg, cmap='YlOrBr', latlon=True)
    elif varname == 'ST': m.pcolormesh(x, y, data_avg, cmap='Set3', latlon=True)
    elif varname == 'SI': m.pcolormesh(x, y, data_avg, cmap='RdYlBu_r', vmax=1, latlon=True)
    elif varname == 'NPP': m.pcolormesh(x, y, data_avg, cmap='viridis', vmin=0, vmax=2000, latlon=True)
    elif varname == 'FPC': m.pcolormesh(x, y, data_avg, cmap='viridis', latlon=True)
    else: print("Variable not found")

    # Add map features
    m.drawcoastlines(linewidth = .3)
    m.drawcountries(linewidth = .3)
    m.drawmapboundary() # alternative background: #m.shadedrelief() or m.bluemarble()
    m.drawparallels(np.arange(lat_min, lat_max + 1, 10), labels=[1, 0, 0, 0], linewidth = .3)
    m.drawmeridians(np.arange(lon_min, lon_max + 1, 10), labels=[0, 0, 0, 1], linewidth = .3)

    if varname=='ST':
        tickloc_arr=np.zeros(13)
        for ind in range(0,13):
            tickloc_arr[ind] = (ind + 0.5)-ind/12
        cbar = m.colorbar(extend='neither', ticks=tickloc_arr)
        labels_arr=np.arange(13)
        cbar.set_ticklabels(labels_arr)
    else:    
        cbar = m.colorbar(extend='neither')

    cbar.set_label(cbar_label, fontsize=14)
    plt.savefig('../Figures/' + varname+'_map.png', dpi=500)
    #plt.title(title, fontsize=16)
    plt.show()




# ***
# #### Load the data:

# In[ ]:

import numpy as np
import netCDF4 as nc
from warnings import filterwarnings
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')

# already masked and averaged over the period
# masked out: all water areas, areas without natural vegetation (urban, glacier, ...), and with too much change in tree cover (deforestation and other).

datafolder="../data_Sebastian"
domain="AmazonSB"

Treecover_data = nc.Dataset(datafolder+'/treecover_'+domain+'_masked.nc', mode='r')
Prec_data = nc.Dataset(datafolder+'/Prec_'+domain+'_masked.nc', mode='r')
Tair_data = nc.Dataset(datafolder+'/Tair_'+domain+'_masked.nc', mode='r')
SWdown_data = nc.Dataset(datafolder+'/SWdown_'+domain+'_masked.nc', mode='r')
LWdown_data = nc.Dataset(datafolder+'/LWdown_'+domain+'_masked.nc', mode='r')
SI_data = nc.Dataset(datafolder+'/SI_'+domain+'_masked.nc', mode='r')
SoilDepth_data = nc.Dataset(datafolder+'/soil_depth_AmazonSB_masked.nc', mode='r')
SoilType_data = nc.Dataset(datafolder+'/soil_type_AmazonSB_masked.nc', mode='r')

TC = np.transpose(Treecover_data.variables['treecover'][:], (0, 2, 1)) # this TC is in (time, lon, lat), so need to switch it to (time, lat, lon) like PR has
TC = np.where(TC>0, TC, np.nan)                                          # removing zeros in non land area (probably sea area); this would not remove the treeless tree cover

PR = Prec_data.variables['Prec'][:]*365*24*3600
TM = Tair_data.variables['Tair'][:]-273.15
SW = SWdown_data.variables['SWdown'][:]
LW = LWdown_data.variables['LWdown'][:]
SI = SI_data.variables['Prec'][:]
SD = SoilDepth_data.variables['soild'][:]
ST = SoilType_data.variables['stexture'][:]


# LPJ model output (unmasked):
file_npp = nc.Dataset(datafolder+'/npp_woody_'+domain+'.nc', mode='r')
NPP_var = file_npp.variables['NPP']
NPP = np.array(NPP_var[:])              # Convert NetCDF variable to a NumPy array
NPP[NPP < 0] = np.nan

file_fpc = nc.Dataset('../data_Sebastian/fpc_woody_'+domain+'.nc')
fpc_var = file_fpc.variables['FPC']
FPC=np.array(fpc_var[:])
FPC[FPC < 0] = np.nan




# ***
# #### Data on the map:

# In[ ]:



lat_values = Treecover_data.variables['lat'][:]
lon_values = Treecover_data.variables['lon'][:]

#"soilmap" : [null,
#  1-3 "clay", "silty clay", "sandy clay", 
#  4-9 "clay loam", "silty clay loam", "sandy clay loam", "loam", "silt loam", "sandy loam", 
# 10-12 "silt", "loamy sand", "sand", 
# 13 "rock and ice"],


#plot_map(TC, lat_values, lon_values, 'TC', 'Tree Cover [%]')
#plot_map(PR, lat_values, lon_values, 'PR', 'Precipitation [mm/year]')
plot_map(TM, lat_values, lon_values, 'TM', 'Temperature [°C]')
plot_map(SW, lat_values, lon_values, 'SW', 'Downwelling short-wave radiation [W/m²]')
plot_map(LW, lat_values, lon_values, 'LW', 'Downwelling long-wave radiation [W/m²]')
#plot_map(SD, lat_values, lon_values, 'SD', 'Soil Depth [m]')
#plot_map(ST, lat_values, lon_values, 'ST', 'Soil Type')
#plot_map(SI, lat_values, lon_values, 'SI', 'Precipitation seasonality')



#plot_map(NPP, lat_values, lon_values, 'NPP', 'NPP [gC/m²/yr]')
#plot_map(FPC, lat_values, lon_values, 'FPC', 'FPC')


