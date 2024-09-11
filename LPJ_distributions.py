#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 18:06:50 2024

@author: bathiany
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib import cm

def transform(data, varname):
    data2 = np.reshape(data, (1,np.product(data.shape)))
    data2 = data2[~np.isnan(data2)]
    data2 = data2[data2<10**19]
    data2 = data2[data2>-10**19]

    if varname=='Tair' or varname=='SWdown' or varname=='LWdown':
        data2 = data2[data2>10]

    if varname=='Prec':
        data2=data2*365*24*3600

    #print(data2.shape)
    return data2


def histogram(data, varname, unit):
    plt.figure()
    dpi=500
    fontsize_ticks=16
    fontsize_title=24
    fontsize_labels=16
    fig, ax = plt.subplots(figsize=(10,7))
    plt.hist(data, bins=100)
    ax.set_xticklabels(fontsize=fontsize_ticks)
    ax.set_xlabel(varname + ' ' + unit, fontsize=fontsize_labels)
    ax.set_ylabel("count", fontsize=fontsize_labels)
    if varname=='Prec':
        plt.xlim(0,4000)

    figname='../Figures/histogram_'+domain+'_'+varname+'.png'
    plt.savefig(figname, dpi=dpi)
    plt.show()


def histogram_partial(data, varname, unit, data2):
    plt.figure()
    dpi=500
    fontsize_ticks=16
    fontsize_title=24
    fontsize_labels=16
    fig, ax = plt.subplots(figsize=(10,7))
    plt.hist(data, bins=100)
    #ax.set_xticklabels(fontsize=fontsize_ticks)
    ax.set_xlabel(varname + ' ' + unit, fontsize=fontsize_labels)
    ax.set_ylabel("count", fontsize=fontsize_labels)

    figname='../Figures/histogram_'+domain+'_partial_'+varname+'.png'
    plt.savefig(figname, dpi=dpi)
    plt.show()


def histogram_colored(data, varname, unit, data2, varname2, unit2):
    plt.figure()
    dpi=500
    fontsize_ticks=16
    fontsize_title=24
    fontsize_labels=16
    fig, ax = plt.subplots(figsize=(10,7))

    #plt.hist(data, bins=100)

    Nbins=100
    Ncells=np.size(data)

    counts, bin_edges = np.histogram(data, bins=Nbins)
    # Calculate the bin centers
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # Set vmin and vmax for the colormap normalization
    if varname2=='Prec':
        vmin = 0   # Minimum value for the color scale
        vmax = 3000   # Maximum value for the color scale
    elif varname2=='FPC_woody':
        vmin = 0   # Minimum value for the color scale
        vmax = 1   # Maximum value for the color scale
    print(vmin)
    # Normalize data2 values within vmin and vmax range
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    
    
    # filter and average data2 to color the bins:
    data2_binned=np.zeros(Nbins)
    for bin in range(0,Nbins):
        data2_binned[bin]=np.mean(data2[(data >= bin_edges[bin]) & (data < bin_edges[bin+1])])
    
    #print(data2_binned)
    
    cmap = cm.viridis
    
    # Create colors based on data2 values
    #colors = cmap(data2_binned)
    colors = cmap(norm(data2_binned))
    
    # Plot histogram with colored bins
    plt.bar(bin_centers, counts/Ncells, width=bin_edges[1] - bin_edges[0], color=colors, edgecolor='black')

    if varname2=='Prec':
        cbar=plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))
        # Change the font size of the colorbar label
        cbar.set_label('Precipitation ' + unit2, fontsize=fontsize_labels)  # Adjust the fontsize as needed
    elif varname2=='FPC_woody':
        cbar=plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))
        cbar.set_label('FPC', fontsize=fontsize_labels)  # Adjust the fontsize as needed


    ## Change the font size of the colorbar tick labels
    cbar.ax.tick_params(labelsize=fontsize_ticks)  # Adjust the fontsize as needed


    #ax.set_xticklabels(fontsize=fontsize_ticks)
    ax.set_ylabel("relative frequency", fontsize=fontsize_labels)
    plt.xticks(fontsize=fontsize_ticks)
    plt.yticks(fontsize=fontsize_ticks)

    if varname=='NPP_woody':
        plt.ylim(0,200/Ncells)
        plt.xlim(1,2000)
        ax.set_xlabel('NPP' + ' ' + unit, fontsize=fontsize_labels)

    figname='../Figures/histogram_'+domain+'_colored_'+varname+'_'+varname2+'.png'
    plt.savefig(figname, dpi=dpi)
    plt.show()



## load data from LPJmL 
#exp='sba0080'
#domain='Mesosouthamerica'
domain='AmazonSB'

#varlist=('npp_woody',)  # fpc_woody

#for var in varlist:
#file_npp = nc.Dataset('/home/bathiany/Projects/Vegetation_resilience_indicators/LPJ/'+exp+'/npp_woody_timmean.nc')
#file_fpc = nc.Dataset('/home/bathiany/Projects/Vegetation_resilience_indicators/LPJ/'+exp+'/fpc_woody_timmean.nc')
file_npp = nc.Dataset('../data_Sebastian/npp_woody_'+domain+'.nc')
file_fpc = nc.Dataset('../data_Sebastian/fpc_woody_'+domain+'.nc')

data_npp = file_npp['NPP']
data_fpc = file_fpc['FPC']


file_Tair=nc.Dataset('../data_Sebastian/Tair_'+domain+'.nc')
file_Prec=nc.Dataset('../data_Sebastian/Prec_'+domain+'.nc')
file_SW=nc.Dataset('../data_Sebastian/SWdown_'+domain+'.nc')
file_LW=nc.Dataset('../data_Sebastian/LWdown_'+domain+'.nc')

data_Tair = file_Tair['Tair']
data_Prec = file_Prec['Prec']
data_SW = file_SW['SWdown']
data_LW = file_LW['LWdown']


npp=transform(data_npp, 'npp')
fpc=transform(data_fpc, 'fpc')
LW=transform(data_LW, 'LWdown')
SW=transform(data_SW, 'SWdown')
Prec=transform(data_Prec, 'Prec')
Tair=transform(data_Tair, 'Tair')



#histogram(npp, "NPP_woody","[gC/yr/m²]")
#histogram(fpc, "fpc_woody", "")
#histogram(Tair, "Tair", "[K]")
#histogram(Prec, "Prec", "[mm/yr]")
#histogram(SW, "SWdown", "[W/m²]")
#histogram(LW, "LWdown", "[W/m²]")

histogram_colored(npp, "NPP_woody","[gC/m²/yr]", fpc, 'FPC_woody', '') 
histogram_colored(npp, "NPP_woody","[gC/m²/yr]", Prec, 'Prec', '[mm/yr]') 

