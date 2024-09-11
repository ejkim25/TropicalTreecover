#!/usr/bin/env python
# coding: utf-8

# ***
# #### Load the data:

# In[ ]:

#import scipy.io
import numpy as np
#import pandas as pd
import netCDF4 as nc
import seaborn as sns
import matplotlib.pyplot as plt
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='\`np.bool\` is a deprecated alias')

# already masked and averaged over the period
# masked out: all water areas, areas without natural vegetation (urban, glacier, ...), and with too much change in tree cover (deforestation and other).

datafolder="../data_Sebastian"
domain="AmazonSB"

Treecover_data = nc.Dataset(datafolder+'/treecover_'+domain+'_masked.nc', mode='r')                       # source: MODIS; period: 2000-2020
Prec_data = nc.Dataset(datafolder+'/Prec_WFDE5_CRU_v1.0_sba0080_'+domain+'_masked.nc', mode='r')                             
Tair_data = nc.Dataset(datafolder+'/Tair_WFDE5_CRU_v1.0_sba0080_'+domain+'_masked.nc', mode='r')                             
SWdown_data = nc.Dataset(datafolder+'/SWdown_WFDE5_CRU_v1.0_sba0080_'+domain+'_masked.nc', mode='r')                             
LWdown_data = nc.Dataset(datafolder+'/LWdown_WFDE5_CRU_v1.0_sba0080_'+domain+'_masked.nc', mode='r')                             
SI_data = nc.Dataset(datafolder+'/SI_WFDE5_CRU_v1.0_sba0080_'+domain+'_masked.nc', mode='r')
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



# ***
# #### Data on the map:

# In[ ]:

fig, axes = plt.subplots(4, 2, figsize=(9,12))

sns.histplot(TC.flatten(), kde=True, bins=40, color='black', fill=False, ax=axes[0,0]) # removed zero values
#axes[0,0].tick_params(labelsize=12)
axes[0,0].set_ylim(0,400)
axes[0,0].set_xlabel('Tree Cover [%]', fontsize=10) # weight='bold'

sns.histplot(TM[0,:,:].flatten(), kde=True, bins=40, color='black', fill=False, ax=axes[0,1])
#axes[0,1].set_ylim(0,700)
axes[0,1].set_xlabel('Temperature [°C]', fontsize=10)

sns.histplot(np.where((PR>=0)&(PR<=4000), PR, np.nan).flatten(), kde=True, bins=40, color='black', fill=False, ax=axes[1,0])
#axes[1,0].set_ylim(0,250)
axes[1,0].set_xlabel('Precipitation [mm/year]', fontsize=10)

sns.histplot(SI[0,:,:].flatten(), kde=True, bins=40, color='black', fill=False, ax=axes[1,1])
axes[1,1].set_xlabel('Precipitation seasonality', fontsize=10)

sns.histplot(SW.flatten(), kde=True, bins=40, color='black', fill=False, ax=axes[2,0])
#axes[1,0].set_ylim(0,300)
axes[2,0].set_xlabel('SWdown [W/m²]', fontsize=10)

sns.histplot(LW.flatten(), kde=True, bins=40, color='black', fill=False, ax=axes[2,1])
#axes[2,1].set_ylim(0,400)
axes[2,1].set_xlabel(r'LWdown [W/m²]', fontsize=10, labelpad=1.5)

sns.histplot(SD.flatten(), kde=True, bins=40, color='black', fill=False, ax=axes[3,0])
#axes[3,0].set_ylim(0,700)
axes[3,0].set_xlabel('Soil Depth [m]', fontsize=10)

sns.histplot(ST.flatten(), kde=True, bins=40, color='black', fill=False, ax=axes[3,1])
#axes[3,1].set_ylim(0,850)
axes[3,1].set_xlabel('Soil Type', fontsize=10)


axes[0,0].set_ylabel('Frequency', fontsize=10)
axes[1,0].set_ylabel('Frequency', fontsize=10)
axes[2,0].set_ylabel('Frequency', fontsize=10)
axes[3,0].set_ylabel('Frequency', fontsize=10)
axes[0,1].set_ylabel('Frequency', fontsize=10)
axes[1,1].set_ylabel('Frequency', fontsize=10)
axes[2,1].set_ylabel('Frequency', fontsize=10)
axes[3,1].set_ylabel('Frequency', fontsize=10)

plt.subplots_adjust(hspace=0.3, wspace=0.25, bottom=0)
plt.savefig('../Figures/histograms_4x2.png', dpi=500, bbox_inches = 'tight', pad_inches = 0.1)
plt.show()
