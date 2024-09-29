#!/usr/bin/env python
# coding: utf-8

# ### Case 1 - Outcome:

# In[ ]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(123)

def productivityTemp(temp):  
    betaT = np.where(temp<20, 0.0, np.where(temp>35, 0.0, -0.0175*(temp-20)*(temp-35)))
    return betaT

def productivityPrec(p):                
    alpha = 0.000054           
    c = 138                    
    gamma = 0.00017*9100

    P1 = c*np.exp(gamma/2) # =~300                                       
    P2 = c*np.exp(gamma/2) + np.exp(gamma)/np.sqrt(0.03*alpha) # =~4000 
    
    betaP = np.where(p < P1, 0.0, np.where(p <= P2, 1.0 - 1.0/(1+alpha*((p-P1)/np.exp(gamma))**2), 1.0))
    return betaP

def Teq(Temp, Prec):
    mort = 0.12
    prod = productivityTemp(Temp)*productivityPrec(Prec)
    BM = prod/mort
    T_eq = np.where(BM < 1.0, 0.0, 1.0 - 1.0/BM)
    return T_eq, BM

def PrecTemp_scatter(temp, prec):
    fontsize_labels = 14
    fontsize_ticks = 14
    
    levels = np.linspace(0,1,51)
    CT=plt.contourf(TEMPv, PRECv, TEQ, levels=levels)
    plt.scatter(temp, prec, s=1, color='black')
    plt.xlabel('Temperature [°C]', fontsize=fontsize_labels)
    plt.ylabel('Precipitation [mm/yr]', fontsize=fontsize_labels)
    plt.xlim(20,35)
    plt.ylim(0,3000)
    plt.xticks([20,25,30,35], fontsize=fontsize_ticks)
    plt.yticks([0,1000,2000,3000], fontsize=fontsize_ticks)
    cbar = plt.colorbar(CT, location='right', ticks=[0, 0.25, 0.5, 0.75, 1])
    cbar.set_label('Tree cover fraction', fontsize=fontsize_labels)
    cbar.ax.tick_params(labelsize=fontsize_labels)
     
def TreeCover_dist(treecover):
    fontsize_labels = 14
    fontsize_ticks = 14
    
    sns.histplot(treecover, kde=True, color='black', bins=40, element="bars", alpha=0.5, edgecolor=None)
    plt.xlabel('Tree cover fraction', fontsize=fontsize_labels)
    plt.ylabel('Frequency', fontsize=fontsize_labels)
    plt.xlim(0,1)
    plt.ylim(0,500)
    plt.xticks([0, 0.25, 0.5, 0.75, 1], fontsize=fontsize_ticks)
    plt.yticks(fontsize=fontsize_ticks)

# <Data>    
# Make data for contour
TEMP = np.linspace(20, 35, 100)
PREC = np.linspace(0, 4000, 100)               
TEMPv, PRECv = np.meshgrid(TEMP, PREC)
TEQ, BM = Teq(TEMPv, PRECv)

# Make data for Teq vs productivity (P-T independent)
datasize1 = 3000
meanT1, stdT1 = 27.5, 3
meanP1, stdP1 = 1000, 500  # pick the center where the gradient is sharper
TEMP_beta1 = np.random.normal(meanT1, stdT1, datasize1)   
PREC_beta1 = np.random.normal(meanP1, stdP1, datasize1)   
TEQ_beta1, BM_beta1 = Teq(TEMP_beta1, PREC_beta1)

# Make data for Teq vs productivity (P-T correlated)
datasize2 = 3000
meanT2, stdT2 = 27.5, 3
slope2 = (1500-500)/(35-20)
cst2 = 500-slope2*20        
stdEpsP2 = 200
TEMP_beta2 = np.random.normal(meanT2, stdT2, datasize2)
PREC_beta2 = slope2*TEMP_beta2 + cst2 + np.random.normal(0, stdEpsP2, TEMP_beta2.size)
TEQ_beta2, BM_beta2 = Teq(TEMP_beta2, PREC_beta2)
corr2 = round(np.corrcoef(TEMP_beta2, PREC_beta2)[0,1],2)
print("corr coeff 2:", corr2)
print("a:", round(slope2,2))
print("b:", round(cst2,2))

# <Graph>
plt.figure(figsize=(10,8))

plt.subplot(2,2,1) # a.Color_Independent
PrecTemp_scatter(TEMP_beta1, PREC_beta1) 
plt.text(-0.32, 1.0, 'a', fontsize=22, ha='center', va='center', transform=plt.gca().transAxes) #fontweight='bold'

plt.subplot(2,2,2) # b.Dist_Independent
TreeCover_dist(TEQ_beta1)                
plt.text(-0.22, 1.0, 'b', fontsize=22, ha='center', va='center', transform=plt.gca().transAxes)

plt.subplot(2,2,3) # c.Color_Correlated
PrecTemp_scatter(TEMP_beta2, PREC_beta2) 
plt.text(-0.32, 1.0, 'c', fontsize=22, ha='center', va='center', transform=plt.gca().transAxes)

plt.subplot(2,2,4) # d.Dist_Correlated
TreeCover_dist(TEQ_beta2)                
plt.text(-0.22, 1.0, 'd', fontsize=22, ha='center', va='center', transform=plt.gca().transAxes)

plt.tight_layout(pad=2.0, w_pad=2.0, h_pad=2.0)
plt.savefig('Results_Model2_Case1.png', dpi=500)
plt.show()
# ==========================================================================================================================================


# ***
# ### Case 1 - Sensitivity Analysis:

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from matplotlib.colors import ListedColormap

np.random.seed(123)

def productivityTemp(temp):                                   
    betaT = np.where(temp<20, 0.0, np.where(temp>35, 0.0, -0.0175*(temp-20)*(temp-35)))
    return betaT

def productivityPrec(p):                        
    alpha = 0.000054           
    c = 138                    
    gamma = 0.00017*9100

    P1 = c*np.exp(gamma/2) # =~300                                       
    P2 = c*np.exp(gamma/2) + np.exp(gamma)/np.sqrt(0.03*alpha) # =~4000  
    
    betaP = np.where(p < P1, 0.0, np.where(p <= P2, 1.0 - 1.0/(1+alpha*((p-P1)/np.exp(gamma))**2), 1.0))
    return betaP

def Teq(temp, prec):
    mort = 0.12
    prod = productivityTemp(temp)*productivityPrec(prec)
    BM = prod/mort
    T_eq = np.where(BM < 1.0, 0.0, 1.0 - 1.0/BM)
    return T_eq, BM

# <Graph>
fig, axes = plt.subplots(3, 3, figsize=(15, 15))
plt.subplots_adjust(hspace=0.3, wspace=-0.4)

datasize = 3000
meanT_rng = np.arange(20, 36, 2.5)
meanP_rng = np.arange(1000, 3001, 500)
x_grid = np.linspace(-0.1, 1.1, 1000)

stdT_values = [2, 3, 4]
stdP_values = [300, 500, 1000]

# Variable to track the max nr of modes across all subplots
global_max_modes = 0

# Iterate through the combinations of stdT and stdP
for row, stdT in enumerate(stdT_values):
    for col, stdP in enumerate(stdP_values):
        num_modes = np.zeros((meanT_rng.size, meanP_rng.size), dtype=int)

        # Calculate num_modes for each combination of meanT and meanP
        for i, meanT in enumerate(meanT_rng):
            for j, meanP in enumerate(meanP_rng):
                TEMP_beta = np.random.normal(meanT, stdT, datasize)
                PREC_beta = np.random.normal(meanP, stdP, datasize)  
                TEQ, BM = Teq(TEMP_beta, PREC_beta)
                kde = gaussian_kde(TEQ)
                density_estimation = kde(x_grid)
                peaks, _ = find_peaks(density_estimation, height=0.835) #prominence=0.1, distance=.
                num_modes[i, j] = len(peaks)
        
        # Update the global max number of modes
        global_max_modes = max(global_max_modes, np.max(num_modes))

        # Custom colormap
        if np.max(num_modes) == 1:
            colors = ['#BBF90F']
        elif np.max(num_modes) == 2:
            colors = ['#BBF90F', '#9ACD32']
        elif np.max(num_modes) == 3:
            colors = ['#BBF90F', '#9ACD32', '#008000']
        custom_cmap = ListedColormap(colors)

        # Subplot for each combination of stdT and stdP
        ax = axes[row, col]
        flip_num_modes = num_modes[::-1] # flip
        ax.imshow(flip_num_modes, cmap=custom_cmap)

        nx = flip_num_modes.shape[1]
        ny = flip_num_modes.shape[0]
        ax.vlines(np.arange(nx)+0.5, ymin=-0.5, ymax=ny-0.5, colors='white', linewidths=0.5)
        ax.hlines(np.arange(ny)+0.49, xmin=-0.5, xmax=nx-0.5, colors='white', linewidths=0.5)

        x_labels = ['1000', '1500', '2000', '2500', '3000']
        y_labels = ['35.0', '32.5', '30.0', '27.5', '25.0', '22.5', '20.0']
        ax.set_xticks(np.arange(len(x_labels)))
        ax.set_xticklabels(x_labels, fontsize=10)
        ax.set_yticks(np.arange(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=10)
        ax.set_xlabel('$\mu_P$ [mm/yr]', fontsize=12)
        ax.set_ylabel('$\mu_T$ [°C]', fontsize=12)
        ax.set_title(f'$\sigma_T$ = {stdT}, $\sigma_P$ = {stdP}', fontsize=13)

# legend in the upper-right corner of the entire figure
if global_max_modes == 1:
    colors = ['#BBF90F']
    labels = ['1']
elif global_max_modes == 2:
    colors = ['#BBF90F', '#9ACD32']
    labels = ['1', '2']
elif global_max_modes == 3:
    colors = ['#BBF90F', '#9ACD32', '#008000']
    labels = ['1', '2', '3']

patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, labels)]
fig.legend(handles=patches, title='No. modes', bbox_to_anchor=(0.87, 0.885), loc='upper right', fontsize='small')
plt.savefig('Results_Model2_Case1_Sensitivity.png', dpi=500)
plt.show()
# ==========================================================================================================================================


# ***
# ### Case 2 - Outcome:

# In[ ]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(123)

def productivityTemp(temp):                                  
    betaT = np.where(temp<20, 0.0, np.where(temp>35, 0.0, -0.0175*(temp-20)*(temp-35)))
    return betaT

def productivityTheta(theta):       
    betaTheta = theta                
    return betaTheta

def productivityPrec(p):                        
    alpha = 0.000054           
    c = 138                    
    gamma = 0.00017*9100

    P1 = c*np.exp(gamma/2) # =~300                                       
    P2 = c*np.exp(gamma/2) + np.exp(gamma)/np.sqrt(0.03*alpha) # =~4000  
    
    betaP = np.where(p < P1, 0.0, np.where(p <= P2, 1.0 - 1.0/(1+alpha*((p-P1)/np.exp(gamma))**2), 1.0))
    return betaP

def Teq(temp, theta, prec):
    mort = 0.12
    prod = productivityTemp(temp)*productivityTheta(theta)*productivityPrec(prec)
    BM = prod/mort
    T_eq = np.where(BM < 1.0, 0.0, 1.0 - 1.0/BM)
    return T_eq, BM

def PrecTemp_scatter(temp, prec, TEQv, title):
    fontsize_title = 14
    fontsize_labels = 14
    fontsize_ticks = 14    
    
    levels = np.linspace(0, 1, 51)
    CT = plt.contourf(TEMPv, PRECv, TEQv, levels=levels) 
    plt.scatter(temp, prec, s=1, color='black')
    plt.xlabel('Temperature [°C]', fontsize=fontsize_labels)
    plt.ylabel('Precipitation [mm/yr]', fontsize=fontsize_labels)
    plt.xlim(20,35)
    plt.ylim(0,3000)
    plt.xticks([20, 25, 30, 35], fontsize=fontsize_ticks)
    plt.yticks([0,1000,2000,3000], fontsize=fontsize_ticks)
    plt.title(title, fontsize=fontsize_title)
    cbar = plt.colorbar(CT, location='right', ticks=[0, 0.25, 0.5, 0.75, 1])
    cbar.set_label('Tree cover fraction', fontsize=fontsize_labels)
    cbar.ax.tick_params(labelsize=fontsize_labels)

def PrecTheta_scatter(theta, prec, TEQv, title):
    fontsize_title = 14
    fontsize_labels = 14
    fontsize_ticks = 14
    
    levels = np.linspace(0, 1, 51)
    CT = plt.contourf(TTv2, PRECv2, TEQv, levels=levels) 
    plt.scatter(theta, prec, s=1, color='black')
    plt.xlabel(r'$\theta$', fontsize=fontsize_labels)
    plt.ylabel('Precipitation [mm/yr]', fontsize=fontsize_labels)
    plt.xlim(0,1)
    plt.ylim(0,3000)
    plt.xticks([0, 0.25, 0.5, 0.75, 1], fontsize=fontsize_ticks)
    plt.yticks([0,1000,2000,3000], fontsize=fontsize_ticks)
    plt.title(title, fontsize=fontsize_title)
    cbar = plt.colorbar(CT, location='right', ticks=[0, 0.25, 0.5, 0.75, 1])
    cbar.set_label('Tree cover fraction', fontsize=fontsize_labels)
    cbar.ax.tick_params(labelsize=fontsize_labels)
    
def TreeCover_dist(teq_ltt, teq_htt, title):
    fontsize_title = 14
    fontsize_labels = 14
    fontsize_ticks = 14
    
    sns.histplot(teq_ltt, kde=True, color='black', bins=20, element="bars", alpha=0.5, edgecolor=None, label=r'$\theta_{Low}$')
    sns.histplot(teq_htt, kde=True, color='royalblue', bins=32, element="bars", alpha=0.5, edgecolor=None, label = r'$\theta_{High}$')
    plt.xlabel('Tree cover fraction', fontsize=fontsize_labels)
    plt.ylabel('Frequency', fontsize=fontsize_labels)
    plt.xlim(0,1)
    plt.ylim(0,500)
    plt.xticks([0, 0.25, 0.5, 0.75, 1], fontsize=fontsize_ticks)
    plt.yticks(fontsize=fontsize_ticks)
    plt.title(title, fontsize=fontsize_title)
    plt.legend(bbox_to_anchor=(0.08, 0.95), fontsize=12, frameon=False)

# <Data>
# Make data for contour 1 (X: temperature; Y: precipitation)
TEMP = np.linspace(20, 35, 100)
PREC = np.linspace(0, 4000, 100)               
LTT = np.full((TEMP.size, PREC.size), 0.25) # lower range of theta
HTT = np.full((TEMP.size, PREC.size), 0.75) # higher range of theta
TEMPv, PRECv = np.meshgrid(TEMP, PREC)
TEQ_ltt, BM_ltt = Teq(TEMPv, LTT, PRECv)
TEQ_htt, BM_htt = Teq(TEMPv, HTT, PRECv)

# Make data for contour 2 (X: theta; Y: precipitation)
TT2 = np.linspace(0, 1, 100)
PREC2 = np.linspace(0, 4000, 100)
TEMP2 = np.full((TT2.size, PREC2.size), 28) #---------------------------------> can change the temperature
TTv2, PRECv2 = np.meshgrid(TT2, PREC2)
TEQ2, BM2 = Teq(TEMP2, TTv2, PRECv2)

# 1) Make data for Teq vs productivity (Precip: normal; Temp: normal; Theta: uniform; P-T indep)
datasize1 = 1500
meanT1, stdT1 = 27.5, 3
meanP1, stdP1 = 1500, 500
LTT_beta1 = np.random.uniform(0.2, 0.3, datasize1) # lower range of theta
HTT_beta1 = np.random.uniform(0.7, 0.8, datasize1) # higher range of theta
TT_beta1 = np.concatenate((LTT_beta1, HTT_beta1), axis=0)
TEMP_beta1 = np.random.normal(meanT1, stdT1, datasize1)  
TEMP_beta1_double = np.random.normal(meanT1, stdT1, datasize1*2)
PREC_beta1 = np.random.normal(meanP1, stdP1, datasize1) 
PREC_beta1_double = np.random.normal(meanP1, stdP1, datasize1*2)
TEQ_LTT1, BM_LTT1 = Teq(TEMP_beta1, LTT_beta1, PREC_beta1)
TEQ_HTT1, BM_HTT1 = Teq(TEMP_beta1, HTT_beta1, PREC_beta1)

data1 = np.concatenate((TEQ_LTT1, TEQ_HTT1), axis=None)
kde1 = gaussian_kde(data1)
x_grid1 = np.linspace(-0.1, 1.1, 1000)
density_estimation1 = kde1(x_grid1)
peaks1, _ = find_peaks(density_estimation1, height=0.835)
num_modes1 = len(peaks1)
print("Number of modes for case 1):", num_modes1)

# 2) Make random data for Teq vs productivity (Precip: normal; Temp: normal; Theta: uniform; P-T corr)
datasize2 = 1500
meanT2, stdT2 = 27.5, 3
stdEpsP2 = 300
LTT_beta2 = np.random.uniform(0.2, 0.3, datasize2) # lower range of theta
HTT_beta2 = np.random.uniform(0.7, 0.8, datasize2) # higher range of theta
TT_beta2 = np.concatenate((LTT_beta2, HTT_beta2), axis=0)
TEMP_beta2 = np.random.normal(meanT2, stdT2, datasize2)  
TEMP_beta2_double = np.random.normal(meanT2, stdT2, datasize2*2)
slope2 = (2200-500)/(35-20) # slope btw temp and prec
cst2 = 500-slope2*20
PREC_beta2 = slope2*TEMP_beta2 + cst2 + np.random.normal(0, stdEpsP2, TEMP_beta2.size)
PREC_beta2_double = slope2*TEMP_beta2_double + cst2 + np.random.normal(0, stdEpsP2, TEMP_beta2_double.size)
TEQ_LTT2, BM_LTT2 = Teq(TEMP_beta2, LTT_beta2, PREC_beta2)
TEQ_HTT2, BM_HTT2 = Teq(TEMP_beta2, HTT_beta2, PREC_beta2)

data2 = np.concatenate((TEQ_LTT2, TEQ_HTT2), axis=None)
kde2 = gaussian_kde(data2)
x_grid2 = np.linspace(-0.1, 1.1, 1000)
density_estimation2 = kde2(x_grid2)
peaks2, _ = find_peaks(density_estimation2, height=0.835)
num_modes2 = len(peaks2)
print("Number of modes for case 2):", num_modes2)

corr2 = round(np.corrcoef(TEMP_beta2, PREC_beta2)[0,1],2)
print("corr coeff 2:", corr2)
print("a:", round(slope2,2))
print("b:", round(cst2,2))

# <Graph>
plt.figure(figsize=(10,8))

plt.subplot(2,2,1) # a.Color_Independent
PrecTheta_scatter(TT_beta1, PREC_beta1_double, TEQ2, 'P - T independent') 
plt.text(-0.32, 1.0, 'a', fontsize=22, ha='center', va='center', transform=plt.gca().transAxes) #fontweight='bold'

plt.subplot(2,2,2) # b.Dist_Independent
TreeCover_dist(TEQ_LTT1, TEQ_HTT1, 'P - T independent')                
plt.text(-0.22, 1.0, 'b', fontsize=22, ha='center', va='center', transform=plt.gca().transAxes)

plt.subplot(2,2,3) # c.Color_Correlated
PrecTheta_scatter(TT_beta2, PREC_beta2_double, TEQ2, 'P - T correlated') 
plt.text(-0.32, 1.0, 'c', fontsize=22, ha='center', va='center', transform=plt.gca().transAxes)

plt.subplot(2,2,4) # d.Dist_Correlated
TreeCover_dist(TEQ_LTT2, TEQ_HTT2, 'P - T correlated')               
plt.text(-0.22, 1.0, 'd', fontsize=22, ha='center', va='center', transform=plt.gca().transAxes)

plt.tight_layout(pad=2.0, w_pad=2.0, h_pad=2.0)
plt.savefig('Results_Model2_Case2.png', dpi=500)
plt.show()

# <Graph - supplementary>
plt.figure(figsize=(10,8))

plt.subplot(2,2,1)
PrecTemp_scatter(TEMP_beta1, PREC_beta1, TEQ_ltt, r'$\theta_L$')
plt.subplot(2,2,2)
PrecTemp_scatter(TEMP_beta2, PREC_beta2, TEQ_ltt, r'$\theta_L$')
plt.subplot(2,2,3)
PrecTemp_scatter(TEMP_beta1, PREC_beta1, TEQ_htt, r'$\theta_H$')
plt.subplot(2,2,4)
PrecTemp_scatter(TEMP_beta2, PREC_beta2, TEQ_htt, r'$\theta_H$')

plt.tight_layout(pad=2.0, w_pad=2.0, h_pad=2.0)
plt.show()
# ==========================================================================================================================================


# ***
# ### Case 2 - Sensitivity Analysis:

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from matplotlib.colors import ListedColormap

np.random.seed(123)

def productivityTemp(temp):                                   
    betaT = np.where(temp<20, 0.0, np.where(temp>35, 0.0, -0.0175*(temp-20)*(temp-35)))
    return betaT
    
def productivityTheta(theta):       
    betaTheta = theta                
    return betaTheta

def productivityPrec(p):                        
    alpha = 0.000054           
    c = 138                    
    gamma = 0.00017*9100

    P1 = c*np.exp(gamma/2) # =~300                                       
    P2 = c*np.exp(gamma/2) + np.exp(gamma)/np.sqrt(0.03*alpha) # =~4000  
    
    betaP = np.where(p < P1, 0.0, np.where(p <= P2, 1.0 - 1.0/(1+alpha*((p-P1)/np.exp(gamma))**2), 1.0))
    return betaP

def Teq(temp, theta, prec):
    mort = 0.12
    prod = productivityTemp(temp)*productivityTheta(theta)*productivityPrec(prec)
    BM = prod/mort
    T_eq = np.where(BM < 1.0, 0.0, 1.0 - 1.0/BM)
    return T_eq, BM

# <Graph>
fig, axes = plt.subplots(3, 3, figsize=(15, 15))
plt.subplots_adjust(hspace=0.3, wspace=-0.4)

datasize = 1500
meanT_rng = np.arange(20, 36, 2.5)        
meanP_rng = np.arange(1000, 3001, 500)
x_grid = np.linspace(-0.1, 1.1, 1000)

stdT_values = [3, 5, 7]         # <--------------------------------------------------------------------------------change
stdP_values = [300, 500, 1000]  # <--------------------------------------------------------------------------------change

# Variable to track the max nr of modes across all subplots
global_max_modes = 0

LTT_beta = np.random.uniform(0.2, 0.3, datasize) # lower range of theta
HTT_beta = np.random.uniform(0.7, 0.8, datasize) # higher range of theta

# Iterate through the combinations of stdT and stdP
for row, stdT in enumerate(stdT_values):
    for col, stdP in enumerate(stdP_values):
        num_modes = np.zeros((meanT_rng.size, meanP_rng.size), dtype=int)

        # Calculate num_modes for each combination of meanT and meanP
        for i, meanT in enumerate(meanT_rng):
            for j, meanP in enumerate(meanP_rng):
                TEMP_beta = np.random.normal(meanT, stdT, datasize)  
                PREC_beta = np.random.normal(meanP, stdP, datasize)  
                TEQ_LTT, BM_LTT = Teq(TEMP_beta, LTT_beta, PREC_beta)
                TEQ_HTT, BM_HTT = Teq(TEMP_beta, HTT_beta, PREC_beta)
                TEQ_TT = np.concatenate((TEQ_LTT, TEQ_HTT), axis=None)
                kde = gaussian_kde(TEQ_TT)
                density_estimation = kde(x_grid)
                peaks, _ = find_peaks(density_estimation, height=0.835)  #prominence=0.1, distance=.
                num_modes[i,j] = len(peaks)
        
        # Update the global max nr of modes
        global_max_modes = max(global_max_modes, np.max(num_modes))
    
        # Custom colormap
        if np.max(num_modes) == 1:                     # for 1 mode
            colors = ['#BBF90F']                    
        elif np.max(num_modes) == 2:                   # for 2 modes
            colors = ['#BBF90F', '#9ACD32']            
        elif np.max(num_modes) == 3:                   # for 3 modes
            colors = ['#BBF90F', '#9ACD32', '#008000']
        custom_cmap = ListedColormap(colors)

        # Subplot for each combination of stdT and stdP
        ax = axes[row, col]
        flip_num_modes = num_modes[::-1] # flip
        ax.imshow(flip_num_modes, cmap=custom_cmap)

        nx = flip_num_modes.shape[1]
        ny = flip_num_modes.shape[0]
        ax.vlines(np.arange(nx)+0.5, ymin=-0.5, ymax=ny-0.5, colors='white', linewidths=0.5)
        ax.hlines(np.arange(ny)+0.49, xmin=-0.5, xmax=nx-0.5, colors='white', linewidths=0.5)
        
        x_labels = ['1000', '1500', '2000', '2500', '3000']
        y_labels = ['35.0', '32.5', '30.0', '27.5', '25.0', '22.5', '20.0']
        ax.set_xticks(np.arange(len(x_labels)))
        ax.set_xticklabels(x_labels, fontsize=10)
        ax.set_yticks(np.arange(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=10)
        ax.set_xlabel('$\mu_P$ [mm/yr]', fontsize=12)
        ax.set_ylabel('$\mu_T$ [°C]', fontsize=12)
        ax.set_title(f'$\sigma_T$ = {stdT}, $\sigma_P$ = {stdP}', fontsize=13)
        
# legend in the upper-right corner of the entire figure        
if global_max_modes == 1:
    colors = ['#BBF90F']
    labels = ['1']
elif global_max_modes == 2:
    colors = ['#BBF90F', '#9ACD32']
    labels = ['1', '2']
elif global_max_modes == 3:
    colors = ['#BBF90F', '#9ACD32', '#008000']
    labels = ['1', '2', '3']

patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, labels)]
fig.legend(handles=patches, title='No. modes', bbox_to_anchor=(0.87, 0.885), loc='upper right', fontsize='small') #borderaxespad=0
plt.savefig('Results_Model2_Case2_Sensitivity.png', dpi=500)
plt.show()

