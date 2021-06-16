'''
@autor Ábram Pérez Herrero
@date 28/05/21

Code for plot dR/dM_{1}dM_{2} in function of M_{1} and M_{2}
'''

import seaborn as sns;
import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt
import MassFunction4 as MF
import matplotlib as mpl

#characteristics of the plots.

mpl.rcParams.update({'font.size': 14,'font.family':'serif'})
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 3.5
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 3.5
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)

#First, we estimate the normalization of each mass function

normlog=MF.normalization("LogNormal",sigma=0.6,Mc=30)
normpower=MF.normalization("powerLaw",alpha=2.2,delta=30)
normcc=MF.normalization("criticalCollapse",beta=2.2,Mc=20)

#The three mass function considered, where we put the normalization and the 
# value of the parameters used
def powerLaw(M):
    return MF.powerLaw(M,normpower,delta=30)
def logNormal(M):
    return MF.logNormal(M,normlog,Mc=30)
def criticalCollapse(M):
    return MF.criticalCollapse(M,normcc,Mc=20)

# Merge Rate obtained in the work, equation 3.58
def mergeRate(Mi,Mj,f,function):
    if Mj!=Mi:
        return (3.7e6)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*min(function(Mi)/Mi,function(Mj)/Mj)*(function(Mi)/Mi+function(Mj)/Mj)
    else :
      return (3.7e6)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*(function(Mj)/(2*Mj))*(function(Mj)/(Mj))

# Define the value of the primordial black hole abundance
f=0.1;
# We use a for loop to obtain the merge rate for each value of Mi(k) and Mj l
# and the vector A1-power, A2-log, A3-CC
l=0;
k=0;
A1=np.zeros((10, 10));
A2=np.zeros((10, 10));
A3=np.zeros((10, 10));
for Mi in [6,8,11,15,19,26,34,45,60,80]:
    for Mj in [6,8,11,15,19,26,34,45,60,80]:
        A1[k,l]=mergeRate(Mi,Mj,f,powerLaw)
        A2[k,l]=mergeRate(Mi,Mj,f,logNormal)       
        A3[k,l]=mergeRate(Mi,Mj,f,criticalCollapse)         
        l=l+1
    l=0;
    k=k+1;
    
# Plot section
index=[6,8,11,15,19,26,34,45,60,80]
columns=[6,8,11,15,19,26,34,45,60,80]

f,ax = plt.subplots(1,3, gridspec_kw={'width_ratios':[1,1,1]},sharey=True)
sns.set(font_scale=1)
df=DataFrame(A1, index, columns)
mask = np.zeros_like(df)
mask[np.tril_indices_from(mask,-1)] = True
b1=sns.heatmap(df,cmap="plasma",mask=mask, square=True,cbar=True,ax=ax[0],cbar_kws={"shrink": 0.5})
b1.set_ylabel('$M_{1}/M_{\odot}$',size=9)
b1.set_xlabel('$M_{2}/M_{\odot}$',size=9)
ax[0].invert_yaxis()
ax[0].text(8,12,'Gpc$^{-3}\:$yr$^{-1}\:M_{\\odot}^{-2}$', fontsize=6)

df=DataFrame(A2, index, columns)
mask = np.zeros_like(df)
mask[np.tril_indices_from(mask,-1)] = True
b2=sns.heatmap(df,cmap="plasma",mask=mask, square=True,cbar=True,ax=ax[1],cbar_kws={"shrink": 0.5})
b2.set_ylabel('')
b2.set_xlabel('$M_{2}/M_{\odot}$',size=9)
b2.set_yticks([])
ax[1].invert_yaxis()
ax[1].text(8,12,'Gpc$^{-3}\:$yr$^{-1}\:M_{\\odot}^{-2}$', fontsize=6)
df=DataFrame(A3, index, columns)
mask = np.zeros_like(df)
mask[np.tril_indices_from(mask,-1)] = True
b3=sns.heatmap(df,cmap="plasma",mask=mask, square=True,cbar=True,ax=ax[2],cbar_kws={"label":"$ d\mathrm{R}\:/\:dM_{1}dM_{2}$","shrink": 0.5})
b3.set_ylabel('')
b3.set_xlabel('$M_{2}/M_{\odot}$',size=9)
b3.set_yticks([])
ax[2].invert_yaxis()
ax[2].text(8,12,'Gpc$^{-3}\:$yr$^{-1}\:M_{\\odot}^{-2}$', fontsize=6)


ax[0].set_title("Power-Law",size=12,style='italic')
ax[1].set_title("Log-Normal",size=12,style='italic')
ax[2].set_title("C.Collapse",size=12,style='italic')


#Save the fifure in eps and pdf format.
plt.savefig("Plots/dRdmdmPlot2D.pdf", bbox_inches="tight")
