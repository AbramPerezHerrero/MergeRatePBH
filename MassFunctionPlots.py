# -*- coding: utf-8 -*-
"""
Created on Sun May 23 20:50:10 2021
Graph showing six types of mass functions considered in different PBH models.
It is created in the Plots folder, with the name of MassFunctionPlot.pdf
@author: ABRAM PÉREZ HERRERO
"""
import matplotlib.pyplot as plt
import MassFunctions as MF
import numpy as np
f=np.linspace(1,100,1000)
k=0;
a=np.zeros((1000,1)); 
b=np.zeros((1000,1));  
c=np.zeros((1000,1)); 
d=np.zeros((1000,1)); 
e=np.zeros((1000,1)); 
g=np.zeros((1000,1)); 
for M in f:
    a[k]=MF.powerLaw(M)
    b[k]=MF.logNormal(M,Mc=25,sigma=0.2)
    c[k]=MF.criticalCollapse(M,Mc=30)
    e[k]=MF.powerLawBroken(M,alpha2=3.5,b=0.6)
    d[k]=MF.powerLawPeak(M,peak=0.4)
    g[k]=MF.multipeak(M,peak1=0.2,peak2=0.3,sigma1=2,sigma2=3)
    k=k+1;

fig, axs = plt.subplots(2,3)
axs[0,0].plot(f,a,'b')
axs[0,1].plot(f,b,'k')
axs[0,2].plot(f,c,'r')
axs[1,0].plot(f,e,'g')
axs[1,1].plot(f,d,'c')
axs[1,2].plot(f,g,'m')
axs[0,0].set_title("Power-Law",position=(0.65, 0.75),size=10,style='italic')
axs[0,1].set_title("Log-Normal",position=(0.65, 0.75),size=10,style='italic')
axs[0,2].set_title("Critical\n Collpase",position=(0.65, 0.7),size=10,style='italic')
axs[1,0].set_title("Power-Law \n Broken",position=(0.65, 0.7),size=10,style='italic')
axs[1,1].set_title("Power-Law \n Peak",position=(0.65, 0.7),size=10,style='italic')
axs[1,2].set_title("Multi-Peak",position=(0.65, 0.75),size=10,style='italic')
axs[0,0].set_ylabel("$\psi(M)$")
axs[1,0].set_ylabel("$\psi(M)$")
axs[0,0].set_xlabel("$M/M_{\odot}$")
axs[0,1].set_xlabel("$M/M_{\odot}$")
axs[0,2].set_xlabel("$M/M_{\odot}$")
axs[1,0].set_xlabel("$M/M_{\odot}$")
axs[1,1].set_xlabel("$M/M_{\odot}$")
axs[1,2].set_xlabel("$M/M_{\odot}$")
axs[0,0].get_xaxis().set_ticks([])
axs[0,0].get_yaxis().set_ticks([])
axs[0,1].get_xaxis().set_ticks([])
axs[0,1].get_yaxis().set_ticks([])
axs[0,2].get_xaxis().set_ticks([])
axs[0,2].get_yaxis().set_ticks([])
axs[1,0].get_xaxis().set_ticks([])
axs[1,0].get_yaxis().set_ticks([])
axs[1,1].get_yaxis().set_ticks([])
axs[1,1].get_xaxis().set_ticks([])
axs[1,2].get_yaxis().set_ticks([])
axs[1,2].get_xaxis().set_ticks([])
fig.savefig("MassFunctionPlot.eps")