# -*- coding: utf-8 -*-
"""
Created on 26/05/21

@author: ABRAM PÉREZ HERRERO
"""

import numpy as np
import matplotlib.pyplot as plt


def mergeRateMonochromatic(f,m,sigmaeq):
    '''
    MERGE RATE FOR A PBH WITH A MONOCHROMATIC MASS DISTRIBUTION.

    Parameters
    ----------
    f : FLOAT
       FRACTION OF DARK MATTER THAT PBH CONSTITUTE
    m : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    sigmaeq : FLOAT
        DENSITY PERTURBATION PROVIDED BY THE REST OF DARK MATTER IN THE
        FORMATION. url:http://dx.doi.org/10.1088/1742-6596/1051/1/012010 

    Returns
    -------
    FLOAT
        MERGE RATE FOR PBH BINARIES PER UNIT OF VOLUME AT PRESENT (z=0). IN
        year^{-1} Gpc^{-1}.

    '''
    return (3.69e6)*(f*0.85)**2*(m)**(-32/37)*((f*0.85)**2+sigmaeq**2)**(-21/74)


sigmaeq=0.005;
f=np.linspace(1e-4,1,10000)
plt.plot(f,mergeRateMonochromatic(f,1,sigmaeq),'r',label="M=30 $M_{\odot}$")
plt.plot(f,mergeRateMonochromatic(f,30,sigmaeq),'b',label="M=100 $M_{\odot}$")
plt.plot(f,mergeRateMonochromatic(f,1000,sigmaeq),'k',label="M=1000 $M_{\odot}$")
plt.text(0.1, 1e4, "M=30 $M_{\odot}$", fontsize=10, rotation = 24,color="blue")
plt.text(0.07, 13e4, "M=1 $M_{\odot}$", fontsize=10, rotation = 24,color="red")
plt.text(0.16, 1e3, "M=1000 $M_{\odot}$", fontsize=10, rotation = 24,color="black")
plt.yscale("log")
plt.xscale("log")
plt.xlabel("$f_{PBH}=\Omega_{PBH}\: / \:\Omega_{DM}$")
plt.ylabel("Merge rate [$Gpc^{-3} yr^{-1}$]")
x=np.linspace(0,1,10000)
#If you want to add the merge rate interval given by LIGO/Virgo.
y1=213
plt.fill_between(x, 2, y1,alpha=1,color='mistyrose')
plt.savefig("Plots/MonochromaticMargeRatePlot.pdf")

