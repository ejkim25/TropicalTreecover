#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 18:49:41 2023

@author: bathiany
"""

import numpy as np

def set_fixed_parameters(PFT):
    global allom1, allom2, allom3 
    global fL, fR, fS
    global vscal,  lmro_ratio, lmro_offset
    global reinickerp, k_latosa, wooddens
    global k_mort, mort_max, height_max, fpcmax
    global MAX_ITER, epsilon, CDEBT_MAXLOAN_DEFICIT, CDEBT_MAXLOAN_MASS, NSEG
    global EPSILON_BISECT_x, EPSILON_BISECT_y, MAXITER_BISECT
    global lightextcoeff, sla, CAmax
    global Lsapl, Rsapl, Ssapl, Hsapl, Csapl
    global k_est
    
    
    ## universal model params
    allom1=100.0     # no unit
    allom2=40.0      # no unit
    allom3=0.67      # /*0.5*/   no unit
    reinickerp=1.6   # no unit   
    k_latosa=6e3     # no unit  , from par/pft.js
    wooddens=2e5     # gC/m3    , from include/tree.h
    k_mort=0.2       # no unit
    mort_max=0.03       # 1/yr max mortality rate
    height_max=100.0    # maximum height of trees in m
    
    fpcmax=0.95   # m2/m2  if set to 0.9, everything already goes crazy
    
    MAX_ITER=10
    epsilon=1.0E-7
    
    lmro_ratio=1.0    ## this is the same for all trees! leaf mass to root mass, from par/pft.js
    lmro_offset=0.5   # from par/pft.js
    
    CDEBT_MAXLOAN_DEFICIT=0.8 # maximum loan as a fraction of deficit  # in allocation_tree.c 
    CDEBT_MAXLOAN_MASS=0.2   # maximum loan as a fraction of (sapwood-cdebt)  # in allocation_tree.c 
    NSEG=20 # number of segments (parameter in numerical methods)
    #CDEBT_PAYBACK_RATE=0.2     # in turnover_tree.c 
    
    ### parameters for leftmostzero
    ## default from LPJ:
    EPSILON_BISECT_x=0.001 #accuracy in x
    EPSILON_BISECT_y=1.0e-10 #accuracy in y           
    MAXITER_BISECT=40 # max iter

    vscal=1  # no nitrogen limitation
    
    k_est=0.12
    
    if PFT==1:
        ## TrBE
        fL=0.5          # 1/2 in 1/yr
        fR=0.5          # 1/2 in 1/yr
        fS=0.03333333   # 1/30 in 1/yr
        lai_sapl=1.5    # m2/m2
        wood_sapl=1.2   #
        #longevity=1.6     
        CAmax=25.0      # m2
        lightextcoeff=0.5
        sla=0.01986     # m2/gC => 500g/m2   in m2/m2


    Lsapl=np.power(lai_sapl*allom1*np.power(wood_sapl,reinickerp)*np.power(4.0*sla/np.pi/k_latosa,reinickerp*0.5)/sla,2.0/(2.0-reinickerp));

    stemdiam=wood_sapl*np.sqrt(4.0*Lsapl*sla/np.pi/k_latosa)
    height_sapl=allom2*np.power(stemdiam,allom3)

    Ssapl=wooddens*height_sapl*Lsapl*sla/k_latosa
    Hsapl=(wood_sapl-1.0)*Ssapl
    Rsapl=(1.0/lmro_ratio)*Lsapl
    Csapl=Lsapl+Rsapl+Ssapl+Hsapl

