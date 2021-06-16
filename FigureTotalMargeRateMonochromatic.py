# -*- coding: utf-8 -*-
"""
Code that creates the figure 3.2 of the work, i.e. 
the collapse rate with monochromatic mass function.

@author: Abram Pérez Herrero 
@Date:14/06/21
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#Plot options
mpl.rcParams.update({'font.size': 14,'font.family':'serif'})

mpl.rcParams['font.family'] = 'serif'

mpl.rc('text', usetex=True)

def mergeRateMono(f,m):
    '''
    Merge rate given by the equation 3.42

    Parameters
    ----------
    f : Float
       primordial black hole abundance. It is the fraction 
       \Omega_{PBH}/\Omega_{DM}
    M: Float
        Primordial black hole mass in solar units.

    Returns
    -------
    FLOAT
      The merge rate of PBH with monochromatic masss function.      

    '''    
    return (3.69e6)*(f*0.85)**2*(m)**(-32/37)*((f*0.85)**2+0.005**2)**(-21/74)
#Plot save and create
f=np.linspace(1e-4,1,10000)
plt.plot(f,mergeRateMono(f,1),'r',label="M=30 $M_{\odot}$")
plt.plot(f,mergeRateMono(f,30),'b',label="M=100 $M_{\odot}$")
plt.plot(f,mergeRateMono(f,1000),'k',label="M=1000 $M_{\odot}$")
plt.text(0.1, 1e4, "M=30 $M_{\odot}$", fontsize=10, rotation = 24,color="blue")
plt.text(0.07, 13e4, "M=1 $M_{\odot}$", fontsize=10, rotation = 24,color="red")
plt.text(0.16, 1e3, "M=1000 $M_{\odot}$", fontsize=10, rotation = 24,color="black")
plt.yscale("log")
plt.xscale("log")
plt.xlabel("$f_{PBH}=\Omega_{PBH}\: / \:\Omega_{DM}$")
plt.ylabel("Merge rate [$Gpc^{-3} yr^{-1}$]")
x=np.linspace(0,1,10000)

#Optional merge rate interval given by LIGO/Virgo
#y1=300
#plt.fill_between(x, 0, y1,alpha=1,color='mistyrose')

plt.savefig("Plots/MergerateMonocrhomatic.pdf")

