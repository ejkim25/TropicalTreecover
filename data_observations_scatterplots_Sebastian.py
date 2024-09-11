#!/usr/bin/env python
# coding: utf-8

# ***
# #### Load the data:

import seaborn as sns
import numpy as np
#import pandas as pd
import netCDF4 as nc
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

SD3 = SoilDepth_data.variables['soild'][:]
SD=np.zeros((np.shape(PR)))
SD[0,:,:]=SD3
SD[SD<0]=float("nan")

ST3 = SoilType_data.variables['stexture'][:]
ST=np.zeros((np.shape(PR)))
ST[0,:,:]=ST3
ST[ST<0]=float("nan")


def figure_pair(data, varname_short, varname, unit):
    plt.figure(figsize=(8,3))
    plt.subplot(1,2,1)
    #print(np.shape(data[(PR>=Pmin)&(PR<=Pmax)]))
    sns.histplot(data[(PR>=Pmin)&(PR<=Pmax)].flatten(), kde=True, color='black', bins=40, fill=False)
    #sns.histplot(data[(PR>=Pmin)&(PR<=Pmax)], kde=True, color='black', bins=40, fill=False)

    #plt.title('distribution for', fontsize=11)
    plt.xlabel(varname + ' ' + unit, weight='bold')
    plt.ylabel('Frequency', weight='bold')
    #plt.ylim(0,60)
    plt.subplot(1,2,2)
    plt.scatter(data[(PR>=Pmin)&(PR<=Pmax)], TC[(PR>=Pmin)&(PR<=Pmax)], s=1, color='black')
    #plt.title('TC vs. T for Pmin≤P≤Pmaxmm/yr', fontsize=11)
    plt.ylim(-5,105)
    plt.xlabel(varname + ' ' + unit, weight='bold')
    plt.ylabel('Tree Cover [%]', weight='bold')
    plt.tight_layout()
    plt.savefig('../Figures/pdf_and_scatter_Pslice_'+varname_short+'.png', dpi=500, bbox_inches = 'tight', pad_inches = 0.1)
    plt.show()


# ***
# In[ ]:
Pmin=1600
Pmax=2200

#Pmin=1700
#Pmax=2100
#
#Pmin=0
#Pmax=4000

#TC_selected = TC[(PR>=Pmin)&(PR<=Pmax)]
#PR_selected = PR[(PR>=Pmin)&(PR<=Pmax)]
#TM_selected = TM[(PR>=Pmin)&(PR<=Pmax)]
##SD_selected = SD[(PR[0,:,:]>=Pmin)&(PR[0,:,:]<=Pmax)]
#
#TC_flattened = TC_selected.flatten()
#PR_flattened = PR_selected.flatten()
#TM_flattened = TM_selected.flatten()
##SD_flattened = SD_selected.flatten()
#TC_flattened_original = TC[(TC>=0)&(TC<=100)].flatten()
#PR_flattened_original = PR[(TC>=0)&(TC<=100)].flatten()

#print(np.shape(TM_selected))
#print(np.shape(TM_flattened))

#print(np.shape(SD))

#print(np.shape(TM.flattened))
#testTM=TM[(PR>=Pmin)&(PR<=Pmax)]
#print(np.shape(TM))
#print(np.shape(testTM))


#testSD=SD[(PR>=Pmin)&(PR<=Pmax)]

#print(np.shape(TM[(PR>=Pmin)&(PR<=Pmax)]))
#print(np.shape(SD[(PR>=Pmin)&(PR<=Pmax)]))


#figure_pair(TM,'TM','Temperature', '[°C]')
#figure_pair(SW,'SW','SWdown', '[W/m²]')
#figure_pair(LW,'LW','LWdown', '[W/m²]')
#figure_pair(SI,'SI','Precipitation seasonality', '')
#figure_pair(SD,'SD','Soil depth', 'm')
figure_pair(ST,'ST','Soil type', '')

