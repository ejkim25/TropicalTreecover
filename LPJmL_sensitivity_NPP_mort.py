#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:02:16 2024

@author: bathiany
"""

import matplotlib.pyplot as plt
import LPJ_EunJoo_functions
import numpy as np
import par


def fig_ts(ts, ts_name):
    plt.figure(dpi=dpi)
    plt.rcParams.update({'font.size': 22})
    fig, ax = plt.subplots(figsize=(8,6))
    plt.plot(ts, color="black", linewidth=1.5)

    if ts_name=="Cperm2":
        plt.ylabel("$C \; [g/m^2]$")
    #    filename_AC=exp + '_' + version + 'S'+str(shuffle)+'_'+parameter+'_vec_AC10_Cperm2_diag.npy'
        #data=np.load(filename_AC)
                                            
    elif ts_name=="C":
        plt.ylabel("$C \; [g/ind]$")
    elif ts_name=="N" or ts_name=="N_LPJ":
        plt.ylabel("$N \; [ind/m^2]$")
    elif ts_name=="fpc":
        plt.ylabel("$fpc \; [m^2/m^2]$")
    elif ts_name=="est":
        plt.ylabel("$establishment \; [ind/m^2yr]$")
    elif ts_name=="mort":
        plt.ylabel("$mortality rate \; [1/yr]$")
    else:
        plt.ylabel(ts_name)

    plt.xlabel('$year$')
    plt.title('test', fontsize=18)
    file_fig='test_'+ts_name+'_ts.png'
    plt.savefig(file_fig)     



def fig_heatmap_bm_vs_mort(data,varname):
    
    plt.figure()
    fontsize_ticks=16
    fontsize_title=24
    fontsize_labels=20
    fig, ax = plt.subplots(figsize=(14,10))
    
    if varname=='Cperm2':
        unit='gC/m2'
        vmin=0
        vmax=10000
    elif varname=='fpc':
        unit='foliage projected cover'
        vmin=0
        vmax=1
 
    Panel=plt.imshow(np.transpose(data), cmap="viridis", vmin=vmin, vmax=vmax)
   
    ax.invert_yaxis()
    ax.set_xticks(np.arange(0,len(bm_forcing)))
    ax.set_yticks(np.arange(0,len(mort_forcing)))
    
    ax.set_xticklabels(bm_forcing_int, fontsize=fontsize_ticks)
    ax.set_yticklabels(mort_forcing, fontsize=fontsize_ticks)
    cbar = plt.colorbar(Panel, fraction=0.030, pad=0.04)
    
    
    #plt.title(varname + ', NPP per ind vs mort', fontsize=fontsize_title)
    cbar.set_label(unit, rotation=90, fontsize=fontsize_labels)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
 
    
    if varname=='fpc':
        _cs2 = ax.contour(np.transpose(data), levels=[0.2, 0.4, 0.6, 0.8] ,origin='lower', colors=['white','orange','red','black'])
        cbar.add_lines(_cs2)

    ax.set_xlabel("productivity [gC/yr/m²]", fontsize=fontsize_labels)
    ax.set_ylabel("mortality rate [1/yr]", fontsize=fontsize_labels)
    
    figname='LPJmL_'+varname+'_bm_vs_mort_heatmap.png'
    
    ax.tick_params(axis='x', rotation=45)
    plt.savefig(figname, dpi=dpi)
    plt.show()



dpi=500

PFT=1
par.set_fixed_parameters(PFT)

## time control
YRS=300

## params

## inputs:
bm_min=50
#bm_max=30000
bm_max=1000
Nvalues_bm=20
bm_forcing=np.linspace(bm_min, bm_max, Nvalues_bm)
bm_forcing_int=np.rint(bm_forcing).astype(int)



wscal_forcing=0.8

mort_min=0
mort_max=0.24
#mort_min=-10
#mort_max=-5

#mort_min=0.0151
#mort_max=0.0153

Nvalues_mort=13
mort_forcing=np.linspace(mort_min, mort_max, Nvalues_mort)

#print(mort_forcing)
#print(bm_forcing)


## prepare matrices for parameter-dependent states
fpc_out=np.zeros((Nvalues_bm, Nvalues_mort))
N_out=np.zeros((Nvalues_bm, Nvalues_mort))
Cperm2_out=np.zeros((Nvalues_bm, Nvalues_mort))


## inicond:
L, R, S, H, D, N, height = 0, 0, 0, 0, 0, 0, 0

## prepare output vectors, ts
fpc_yearly=np.zeros(YRS)
N_yearly=np.zeros(YRS)
Cperm2_yearly=np.zeros(YRS)

### loop over parameter index

for cellind_bm in range(0,Nvalues_bm):
    for cellind_mort in range(0,Nvalues_mort):
        
        for year in range(1,YRS):
        
            ## run reduced LPJmL:          
            L, R, S, H, D, N, fpc, height, mort = LPJ_EunJoo_functions.LPJ_Cbalance_A1N1(L, R, S, H, D, N, height, bm_forcing[cellind_bm], wscal_forcing, mort_forcing[cellind_mort])

            #print(L, R, S, H, D, N, fpc, height, mort)
        
            fpc_yearly[year]=fpc
            N_yearly[year]=N
            Cperm2_yearly[year]=N*(L+R+S+H)
        
        #fig_ts(fpc_yearly, 'fpc')
        #fig_ts(N_yearly, 'N')
        #fig_ts(Cperm2_yearly, 'Cperm2')

        fpc_out[cellind_bm, cellind_mort]=fpc
        N_out[cellind_bm, cellind_mort]=N
        Cperm2_out[cellind_bm, cellind_mort]=N*(L+R+S+H)

#fig_heatmap_bm_vs_mort(Cperm2_out,'Cperm2')
fig_heatmap_bm_vs_mort(fpc_out,'fpc')
