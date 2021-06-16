
"""
Created on Sun May 23 20:50:10 2021
Graph showing six types of mass functions considered in different PBH models.
It is created in the Plots folder, with the name of MassFunctionPlot.pdf
@author: ABRAM PÉREZ HERRERO
"""
import matplotlib.pyplot as plt
import MassFunction4 as MF
import numpy as np
import matplotlib as mpl
#Mass interval to plot
Mx=np.linspace(3,100,1000)
k=0 # counter
#Plot options
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

mpl.rcParams['legend.edgecolor'] = 'inherit'
#Arrays to stored the values of each Mass function
a1=np.zeros((1000,1));
a2=np.zeros((1000,1)); 
a3=np.zeros((1000,1)); 
b1=np.zeros((1000,1));  
b2=np.zeros((1000,1));  
b3=np.zeros((1000,1));  
c1=np.zeros((1000,1)); 
c2=np.zeros((1000,1)); 
c3=np.zeros((1000,1)); 
d1=np.zeros((1000,1)); 
for M in Mx:
    a1[k]=MF.powerLaw(M,delta=20,alpha=4,A=MF.normalization("powerLaw",alpha=4,delta=20))
    a2[k]=MF.powerLaw(M,delta=20,alpha=2.2,A=MF.normalization("powerLaw",alpha=2.2,delta=20))
    a3[k]=MF.powerLaw(M,delta=20,alpha=2.7,A=MF.normalization("powerLaw",alpha=2.7,delta=20))
    b1[k]=MF.logNormal(M,Mc=30,sigma=0.4,A=MF.normalization("LogNormal",sigma=0.4,Mc=30))
    b2[k]=MF.logNormal(M,Mc=30,sigma=0.6,A=MF.normalization("LogNormal",sigma=0.6,Mc=30))
    b3[k]=MF.logNormal(M,Mc=30,sigma=1,A=MF.normalization("LogNormal",sigma=1,Mc=30))
    c1[k]=MF.criticalCollapse(M,Mc=30,beta=2.85,A=MF.normalization("criticalCollapse",beta=2.85,Mc=30))
    c2[k]=MF.criticalCollapse(M,Mc=30,beta=1.2,A=MF.normalization("criticalCollapse",beta=1.2,Mc=30))
    c3[k]=MF.criticalCollapse(M,Mc=30,beta=1.7,A=MF.normalization("criticalCollapse",beta=1.7,Mc=30))
    d1[k]=MF.logNormal(M,Mc=25,sigma=0.01,A=MF.normalization("LogNormal",sigma=0.01,Mc=30))
    k=k+1;
    
#Plot save
fig, axs = plt.subplots(1,3)

axs[0].plot(Mx,a1,'k',label="$\\alpha=4$")
axs[0].plot(Mx,a2,'m',label="$\\alpha=2.2$")
axs[0].plot(Mx,a3,'c',label="$\\alpha=2.7$")
axs[1].plot(Mx,b1,'k',label="$\sigma=0.4$")
axs[1].plot(Mx,b2,'c',label="$\sigma=0.6$")
axs[1].plot(Mx,b3,'m',label="$\sigma=1$")
axs[2].plot(Mx,c1,'k',label="$\\beta=2.85$")
axs[2].plot(Mx,c2,'c',label="$\\beta=1.2$")
axs[2].plot(Mx,c3,'m',label="$\\beta=1.7$")

axs[0].set_title("Power-Law",size=10,style='italic')
axs[1].set_title("Log-Normal",size=10,style='italic')
axs[2].set_title("C.Collapse",size=10,style='italic')
axs[0].legend(loc="upper right",prop={'size': 8})
axs[1].legend(loc="upper right",prop={'size': 8})
axs[2].legend(loc="upper right",prop={'size': 8})



axs[0].set_ylabel("$\psi(M)$",size=9)
axs[1].set_ylabel("$\psi(M)$",size=9)
axs[2].set_ylabel("$\psi(M)$",size=9)


axs[0].set_xlabel("$M/M_{\odot}$",size=9)
axs[1].set_xlabel("$M/M_{\odot}$",size=9)
axs[2].set_xlabel("$M/M_{\odot}$",size=9)


axs[0].get_xaxis().set_ticks([])
axs[0].get_yaxis().set_ticks([])
axs[1].get_xaxis().set_ticks([])
axs[1].get_yaxis().set_ticks([])
axs[2].get_xaxis().set_ticks([])
axs[2].get_yaxis().set_ticks([])

plt.subplots_adjust(left=0.125,
                    bottom=0.6, 
                    right=1.2, 
                    top=1, 
                    wspace=0.2, 
                    hspace=0.4)

plt.savefig("Plots/MassFunctionPlotselected.pdf", bbox_inches="tight")
#fig.savefig("Plots/MassFunctionPlot.eps")
# -*- coding: utf-8 -*-


