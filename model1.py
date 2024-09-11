#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def Teq(x1, x2):
    Teq1 = b11/(1.0+np.exp(b12-b13*x1))+b14
    Teq2 = b21/(1.0+np.exp(b22-b23*x2))+b24
    Teqf = Teq1*Teq2  
    return Teqf, Teq1, Teq2
#     
#def X1X2_scatter(ax, x1, x2, title, fontsize=11):
#    levels = np.linspace(0,1,51)
#    #CT=ax.contourf(X1v, X2v, TEQ, 20) # old
#    CT=ax.contourf(X1v, X2v, TEQ, levels=levels)
#    ax.scatter(x1, x2, s=1, color='black')
#    ax.set_xlabel('$x_1$', fontsize=fontsize)
#    ax.set_ylabel('$x_2$', fontsize=fontsize)
#    ax.set_title(title, fontsize=fontsize)
#    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
#    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
#    #ax.set_xlim(0,1)
#    #ax.set_ylim(0,1)
#    ax.set_xlim(-1,1)
#    ax.set_ylim(-1,1)
#    plt.colorbar(CT, location='right', label='Tree Cover fraction', ticks=[0, 0.25, 0.5, 0.75, 1])
#

version=2


## old
#b11, b12, b13, b14 = 0.88, 5, 16, 0.05   # Teq1
#if version==1:
#    b21, b22, b23, b24 = 0.62, 5, 11, 0.35   # Teq2: for distinct sigmoid curves
#elif version==2:
#    b21, b22, b23, b24 = 0.88, 5, 11, 0.05  # Teq2: for close sigmoid curves

## modify to shift left:
b11, b12, b13, b14 = 0.88, -3, 16, 0.05   # Teq1
if version==1:
    b21, b22, b23, b24 = 0.88, -0.5, 11, 0.05   # Teq2: for close sigmoid curves
elif version==2:
    b21, b22, b23, b24 = 0.62, -0.5, 11, 0.35   # Teq2: for distinct sigmoid curves



# Data ===========================================================================    
# Make data for contour
#X1 = np.linspace(0, 1, 100)
#X2 = np.linspace(0, 1, 100)               
X1 = np.linspace(-1, 1, 100)
X2 = np.linspace(-1, 1, 100)   
X1v, X2v = np.meshgrid(X1, X2)
TEQ, TEQ1, TEQ2 = Teq(X1v, X2v)

# Make data for Teq vs variable (normal)
N = 2000
#meanX1, stdX1 = 0.5, 0.3
#meanX2, stdX2 = 0.5, 0.3

meanX1, stdX1 = 0, 0.3
meanX2, stdX2 = 0, 0.3

dataX1 = np.random.normal(meanX1, stdX1, N)   
dataX2 = np.random.normal(meanX2, stdX2, N)
dataTEQ, dataTEQ1, dataTEQ2 = Teq(dataX1, dataX2)


dpi=500

# Graph ==========================================================================
# X1 vs. X2 scatter
#fig, ((ax)) = plt.subplots(figsize=(4,3), nrows=1, ncols=1)
#X1X2_scatter(ax1, dataX1, dataX2, '$X_1, X_2$~$N$({}, {}\u00b2); $X_1 \perp X_2$'.format(meanX1, stdX1))
#X1X2_scatter(ax1, dataX1, dataX2, '')
plt.figure(figsize=(4,3), constrained_layout=True)
fontsize=11
levels = np.linspace(0,1,51)
#CT=ax.contourf(X1v, X2v, TEQ, 20) # old
CT=plt.contourf(X1v, X2v, TEQ, levels=levels)
plt.scatter(dataX1, dataX2, s=1, color='black')
plt.xlabel('$x_1$', fontsize=fontsize)
plt.ylabel('$x_2$', fontsize=fontsize)
plt.xticks([-1, -0.5, 0, 0.5, 1])
plt.yticks([-1, -0.5, 0, 0.5, 1])
#ax.set_xlim(0,1)
#ax.set_ylim(0,1)
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.colorbar(CT, location='right', label='Tree Cover fraction', ticks=[0, 0.25, 0.5, 0.75, 1])
plt.savefig('../Figures/Results_Model1_Sigm_Color'+str(version)+'.png', dpi=dpi)
plt.show()

#Figures/Results_Model1_Sigm_Color2.png

# Sigmoid function
#X1Range=np.linspace(0,1,100)
#X2Range=np.linspace(0,1,100)
X1Range=np.linspace(-1,1,100)
X2Range=np.linspace(-1,1,100)

# curves
TCF_fct, TC1_fct, TC2_fct = Teq(X1Range, X2Range)
plt.figure(figsize=(4,3), constrained_layout=True)
plt.plot(X1Range, TC1_fct, color='blue', label='$\Psi_1$')
plt.plot(X2Range, TC2_fct, color='black', label='$\Psi_2$')
plt.xlabel('$x$', fontsize=12)
plt.ylabel('Tree cover fraction', fontsize=12)
plt.legend(loc='lower right')
plt.ylim(-0.05, 1.05)
plt.yticks([0, 0.25, 0.5, 0.75, 1])
plt.savefig('../Figures/Results_Model1_Sigm_Curve'+str(version)+'.png', dpi=dpi)
plt.show()

# Tree cover distribution w/ kde
fig = plt.figure(figsize=(4,3), constrained_layout=True)
sns.histplot(dataTEQ, kde=True, color='black', bins=20, fill=True, element="step", edgecolor=None)
plt.xlim(0,1)
plt.ylim(0,300)
plt.xlabel('Tree cover fraction', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.xticks([0, 0.25, 0.5, 0.75, 1])
#plt.tight_layout()
plt.savefig('../Figures/Results_Model1_Sigm_Dist'+str(version)+'.png', dpi=dpi)
plt.show()

