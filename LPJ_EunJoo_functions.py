#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 18:36:57 2023

@author: bathiany
"""
import numpy as np
import par


### read LPJ output files are masked arrays,
## convert to normal ones:
def unmask(data):
    data2 = np.float32(data)
    data_unmasked = np.ma.filled(data2, fill_value=0.)
    return data_unmasked
    



###### individual modules (processes)

def fpc_tree(L,N,CAind):
    if CAind>0.0:
        LAIind=lai_tree(L, CAind)
        FPCind=1.0-np.exp(-par.lightextcoeff*LAIind)
        FPC_PFT=CAind*N*FPCind
    else:
        FPC_PFT=0.0
    return FPC_PFT



def lai_tree(L, CAind):
    if CAind>0:
        LAIind=L*par.sla/CAind
    else:
        LAIind=0
    return LAIind



def allometry_tree(L,S,H):

    if ( S<=0.0 or L<=0.0):
        height=0
    else:
        height=S*par.k_latosa/(par.wooddens*L*par.sla)

    if (height>par.height_max):
        height=par.height_max
        sm_ind_temp=S
        S=L*par.height_max*par.wooddens*par.sla/par.k_latosa
        H=H+sm_ind_temp-S

    allometry=par.allom1*np.power(height/par.allom2,par.reinickerp/par.allom3)
    CAind=min(allometry,par.CAmax)
    return S, H, CAind, height


def fcn(leaf_inc,k1,lm,k3,b,L,H):
    return k1*(b-leaf_inc*lm+H)-np.power((b-leaf_inc*lm)/(L+leaf_inc)*k3,1.0+2.0/par.allom3)
    

def bisect(func, xlow, xhigh, k1,lm,k3,b,L,H, xacc, yacc, maxit):
    ymin=1e09     ## see numeric/bisect.c    
    ylow=func(xlow, k1,lm,k3,b,L,H)
    i=0
    while i < maxit:
        xmid=(xlow+xhigh)*0.5
        if xhigh-xlow<xacc:
            return xmid
        
        ymid=func(xmid,k1,lm,k3,b,L,H)
        if(abs(ymid)<ymin):
            ymin=abs(ymid)
            xmin=xmid
        
        if(abs(ymid)<yacc):
            return xmid
        
        if(ylow*ymid<=0):
            xhigh=xmid
        else:
            xlow=xmid
            ylow=ymid
        i=i+1
    
    return xmin

def leftmostzero(func,x1, x2, k1,lm,k3,b,L,H, xacc, yacc, maxiter):
    if(x2<x1):
        swap=x1
        x1=x2
        x2=swap

    dx=(x2-x1)/par.NSEG
    if(func(x1,k1,lm,k3,b,L,H)<0):
         xmid=x1+dx
         while func(xmid,k1,lm,k3,b,L,H)<0 and xmid<=x2-dx:
              xmid=xmid+dx
    else:  #func>0
        xmid=x1+dx
        while func(xmid,k1,lm,k3,b,L,H)>0 and xmid<=x2-dx:
            xmid=xmid+dx
            
    return bisect(func,xmid-dx,xmid,k1,lm,k3,b,L,H,xacc,yacc,maxiter)




#def allocation_tree(L, R, S, H, D, height, bm_inc_ind, wscal):
def allocation_tree(L, R, S, H, D, N, height, bm_inc, wscal):

    #drought1, drought2 = 0, 0
    lmtorm=par.lmro_ratio*(par.lmro_offset+(1-par.lmro_offset)*min(par.vscal,wscal))

    if N<10**-15:
        bm_inc_ind=1
    else:
        bm_inc_ind=bm_inc/N
    #### This was feeding in bm_inc from LPJ and divide by N using the state of N in this reduced model.
#        # not used here anymore. It was useful to make the model stable when feeding in the four 
#        # allocated fluxes (model10) because otherwise becomes numerically unstable.
#                              
#    ## alternative: feed in bm_inc_ind directly, i.e. the NPP per individual,
#    ## then the model can decide how many trees there are, and 
#    ## hence also how large the overall NPP per m2 is.
    
    tinc_H=0
    tinc_D=0
    if lmtorm<1.0e-10:
        print("missing if case")
    # else

    if height>0:
        tinc_ind_min_leaf=par.k_latosa*S/(par.wooddens*height*par.sla)-L
        tinc_ind_min_root=par.k_latosa*S/(par.wooddens*height*par.sla*lmtorm)-R
    else:
        tinc_ind_min_leaf=0
        tinc_ind_min_root=0
        
    cmass_deficit=tinc_ind_min_leaf+tinc_ind_min_root-bm_inc_ind


    if cmass_deficit>0:
        cmass_loan=max(min(cmass_deficit*par.CDEBT_MAXLOAN_DEFICIT,S-D)*par.CDEBT_MAXLOAN_MASS,0)
        bm_inc_ind=bm_inc_ind+cmass_loan
        tinc_D=cmass_loan
        
    else:    # true, except at dry cells!
        tinc_D=0 
    
    
    if (tinc_ind_min_root>=0 and tinc_ind_min_leaf>=0 and (tinc_ind_min_root+tinc_ind_min_leaf<=bm_inc_ind or bm_inc_ind<=0)):
        
        b=S+bm_inc_ind-L/lmtorm+R
        lm=1+1/lmtorm
        k1=np.power(par.allom2,2.0/par.allom3)*4.0*(1.0/np.pi)/par.wooddens
        k3=par.k_latosa/par.wooddens/par.sla
      
        x2=(bm_inc_ind-(L/lmtorm-R))/lm
        if L<1e-10:
            x1=x2/par.NSEG
        else:
            x1=0

        ## bisection:
        if((x1==0 and x2==0) or b-x1*lm<0 or L+x1<=0 or b-x2*lm<0 or L+x2<=0):
            tinc_L=0
        else:
            tinc_L=leftmostzero(fcn,x1,x2,k1,lm,k3,b,L,H,par.EPSILON_BISECT_x,par.EPSILON_BISECT_y,par.MAXITER_BISECT) # params like in LPJ

        if (tinc_L<0):
            tinc_R=0
            
        else:
            tinc_R=(tinc_L+L)/lmtorm-R 

        if(bm_inc_ind>0 and tinc_R+tinc_L>bm_inc_ind):
            tinc_R=bm_inc_ind*tinc_R/(tinc_R+tinc_L)
            tinc_L=bm_inc_ind*tinc_L/(tinc_R+tinc_L)
        
        tinc_S=bm_inc_ind-tinc_L-tinc_R
        

    else: ##   if NOT: (tinc_ind_min_root>=0 and tinc_ind_min_leaf>=0 and (tinc_ind_min_root+tinc_ind_min_leaf<=bm_inc_ind or bm_inc_ind<=0)):

        ## Abnormal allocation:
        tinc_L=(bm_inc_ind+R-L/lmtorm)/(1+1/lmtorm)
        if (tinc_L>0): 
            ## Sitch 2003:
            # "In years of stress, the biomass increment may not allow sufficient
            # allocation to the leaves to fully utilize the current
            # sapwood (given the constraint implied by Eqn. 1: LA = k_lasa * SA). This
            # year's production is then allocated to leaves and roots
            # only, and the excess sapwood mass transferred to the
            # nonliving heartwood pool."
            # follows the "pipe model"
            # "This relationship is based on numerous studies indicating that each unit of leaf area must
            # be supported by a corresponding area of transport tissue (Shinozaki et al., 1964a, b; Kaufmann & Troendle, 1981;
            # Waring et al., 1982; Ryan, 1989; Robichaud & Methven, 1992; Berninger & Nikinmaa, 1994)."
            #drought1=1
            tinc_R=bm_inc_ind-tinc_L
        else:          
            ## Sitch 2003:
            #  "In a year with severe drought
            # there may be insufficient biomass increment to maintain
            # both current sapwood mass and leaf mass. In this case all
            # of the biomass increment is allocated to fine roots and
            # excess sapwood (...) transferred to the heart-
            # wood (...) pool."
            #drought2=1
            tinc_R=bm_inc_ind
            tinc_L=(R+tinc_R)*lmtorm-L
            # tinc_L => litter
        # "excess sapwood" because it is too much for the few leaves,
        # based on empirical relationship between leaf to sapwood cross-sect area
        tinc_S=(tinc_L+L)*par.wooddens*height*par.sla/par.k_latosa-S   
        tinc_H=-tinc_S   # "transferred to the heartwood pool"

    return tinc_L, tinc_R, tinc_S, tinc_H, tinc_D





def turnover_tree(L, R, S, H):
    
    turnover_L=-par.fL*L
    turnover_R=-par.fR*R
    turnover_S=-par.fS*S
    turnover_H=par.fS*S
    turnover=turnover_L+turnover_R+turnover_S
    L=L+turnover_L
    R=R+turnover_R
    S=S+turnover_S
    H=H+turnover_H
    L, R, S, H = crop_pools(L, R, S, H)

    return L, R, S, H, turnover



def tree_mortality(Aperind, N, turnover, L, mort_prescribed):

    if mort_prescribed<0:
        bm_delta=Aperind+turnover
        if(bm_delta<0):
            bm_delta=0        

    if L<=0:
        mort=0
        bm_delta=0
    else:
        if mort_prescribed<0:
            mort=par.mort_max/(1+par.k_mort*bm_delta/(L*par.sla))
        else:
            mort=mort_prescribed
            bm_delta=0
    
    Nmort=-N*mort
    N=N+Nmort
    if N<0:
        N=10**-10
        
    return N, mort



def tree_establishment(fpc, L, R, S, H, D, N):
    
    if (fpc >= par.fpcmax):
        S, H, CAind, height = allometry_tree(L,S,H)
        est=0
    else:
        est=par.k_est*(1-np.exp(-5*(1-fpc)))*(1-fpc)

        Nold=N
        Nnew=N+est
        
        #print(Nold, est, Nnew)
        L=(L*Nold+par.Lsapl*est)/Nnew   
        R=(R*Nold+par.Rsapl*est)/Nnew 
        S=(S*Nold+par.Ssapl*est)/Nnew 
        H=(H*Nold+par.Hsapl*est)/Nnew
        D=D*Nold/Nnew
        
        S, H, CAind, height = allometry_tree(L,S,H)
        N=Nnew
        
    fpc=fpc_tree(L,N,CAind)
        
    return N, L, R, S, H, D, CAind, height, fpc, est





def adjust_tree(fpc, L, N, CAind):
    i = 0
    if (fpc > par.fpcmax):

        while par.fpcmax < fpc and i < par.MAX_ITER:
            frac=par.fpcmax/fpc
            N=N*frac
            fpc=fpc_tree(L,N,CAind)
            i = i + 1
        
    return N, fpc



def crop_pools(L, R, S, H):
    if L<0:
        L=10**-10
    if R<0:
        R=10**-10    
    if S<0:
        S=10**-10    
    if H<0:
        H=10**-10

    return L, R, S, H



# interactive allo, interactive dynveg
#def LPJ_Cbalance_A1N1(L, R, S, H, D, N, height, bm_inc_ind, wscal, mort_prescribed):
def LPJ_Cbalance_A1N1(L, R, S, H, D, N, height, bm_inc, wscal, mort_prescribed):
    
    ####### turnover
    L, R, S, H, turnover = turnover_tree(L, R, S, H)
        
    ###### allocation
    #AL, AR, AS, AH, AD = allocation_tree(L, R, S, H, D, height, bm_inc_ind, wscal)
    AL, AR, AS, AH, AD = allocation_tree(L, R, S, H, D, N, height, bm_inc, wscal)

    
    Aperind=AL+AR+AS+AH
    
    L, R, S, H, D = L+AL, R+AR, S+AS, H+AH, D+AD
    L, R, S, H = crop_pools(L, R, S, H)

    S, H, CAind, height=allometry_tree(L,S,H)

    fpc=fpc_tree(L,N,CAind)

    ##### mortality
    N, mort = tree_mortality(Aperind, N, turnover, L, mort_prescribed)    

    fpc=fpc_tree(L,N,CAind)

    
    ######## establishment  (OTHER PFTS MATTER)  
    N, L, R, S, H, D, CAind, height, fpc, est = tree_establishment(fpc, L, R, S, H, D, N)
    
    ####### adjustment (OTHER PFTS MATTER)
    N, fpc = adjust_tree(fpc, L, N, CAind)


    return L, R, S, H, D, N, fpc, height, mort



