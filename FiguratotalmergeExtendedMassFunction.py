# -*- coding: utf-8 -*-
"""
Code for obtain the plot 3.4 of the work
@author: Abram Pérez Herrero
@Date:14/06/21
"""

from scipy.integrate import dblquad as din
import matplotlib.pyplot as plt
import MassFunction4 as MF
import numpy as np
import matplotlib as mpl
#Normalization of the different mass function used
I=MF.normalization("powerLaw",alpha=2.2,delta=20)
C=MF.normalization("LogNormal",sigma=1,Mc=30)
H=MF.normalization("criticalCollapse",beta=2.85,Mc=30)
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

def lognormal(M,I):
    '''
    Mass function used which takes into account normalisation
    and sigma value.

    Parameters
    ----------
    M : FLOAT
        Primordial black hole mass in solar mass units
    I : FLOAT
        Normalisation of the mass function taking into account the parameters.


    Returns
    -------
    FLOAT
        The probability of the mass given.

    '''
    return MF.logNormal(M,A=I,sigma=1,Mc=30)
def mergeRatePLB(Mi,Mj,f,C):
    zeta=1 #Clustering parameter if it is one implies a random sp.distribution
    '''
         Merge rate given by the equation 3.58

    Parameters
    ----------
    Mi : FLOAT
        Correspond to M1 in the work, are the PBH's mass in solar units.
    Mj : FLOAT
        Correspond to M2 in the work, are the PBH's mass in solar units.
    f : FLOAT
        The abundance of primordial black holes. It coulb be in the interval
        [0,1]
    C : FLOAT
        Normalization of the mass function.
    mu : FLOAT
    Mc of the mass function that has been chosen.

    Returns
    -------
    FLOAT
        The merge rate of PBH with extended masss function. 

    '''
    function=lognormal
    if Mj!=Mi:
            return (3.7e6)*zeta**(16/37)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*min(function(Mi,C)/Mi,function(Mj,C)/Mj)*(function(Mi,C)/Mi+function(Mj,C)/Mj)
    else :
            return (3.7e6)*zeta**(16/37)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*(function(Mj,C)/(2*Mj))*(function(Mj,C)/(Mj))

def powerlaw(M,I):
    '''
    Mass function used which takes into account normalisation
    and alpha value.

    Parameters
    ----------
    M : FLOAT
        Primordial black hole mass in solar mass units
    I : FLOAT
        Normalisation of the mass function taking into account the parameters.

    Returns
    -------
    FLOAT
        The probability of the mass given.

    '''
    return MF.powerLaw(M,A=I,alpha=2.2,delta=20)
def mergeRatePLBlaw(Mi,Mj,f,C):
    
        zeta=1
        '''
         Merge rate given by the equation 3.58
 
         Parameters
         ----------
         Mi : FLOAT
             Correspond to M1 in the work, are the PBH's mass in solar units.
        Mj : FLOAT
            Correspond to M2 in the work, are the PBH's mass in solar units.
            f : FLOAT
            The abundance of primordial black holes. It coulb be in the interval
            [0,1]
            C : FLOAT
            Normalization of the mass function.


        Returns
        -------
        FLOAT
        The merge rate of PBH with extended masss function. 
        
        '''
        function=powerlaw
        if Mj!=Mi:
            return (3.7e6)*zeta**(16/37)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*min(function(Mi,C)/Mi,function(Mj,C)/Mj)*(function(Mi,C)/Mi+function(Mj,C)/Mj)
        else :
            return (3.7e6)*zeta**(16/37)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*(function(Mj,C)/(2*Mj))*(function(Mj,C)/(Mj))

def critical(M,I):
    '''
    Mass function used which takes into account normalisation
    and sigma value.It is created in principle so that the normalisation
    integral does not get into the loop.

    Parameters
    ----------
    M : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    I : FLOAT
        Normalisation of the mass function taking into account the parameters.


    Returns
    -------
    FLOAT
        The probability of the mass given.

    '''
    return MF.criticalCollapse(M,A=I,beta=2.85,Mc=30)
def mergeRatePLcc(Mi,Mj,f,C):
    '''
    Merge rate given by the equation 3.58

    Parameters
    ----------
    Mi : FLOAT
        Correspond to M1 in the work, are the PBH's mass in solar units.
    Mj : FLOAT
        Correspond to M2 in the work, are the PBH's mass in solar units.
    f : FLOAT
        The abundance of primordial black holes. It coulb be in the interval
        [0,1]
    C : FLOAT
        Normalization of the mass function.


    Returns
    -------
    FLOAT
        The merge rate of PBH with extended masss function. 

    '''
    zeta=100
    function=critical
    if Mj!=Mi:
            return (3.7e6)*zeta**(16/37)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*min(function(Mi,C)/Mi,function(Mj,C)/Mj)*(function(Mi,C)/Mi+function(Mj,C)/Mj)
    else :
            return (3.7e6)*zeta**(16/37)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*(function(Mj,C)/(2*Mj))*(function(Mj,C)/(Mj))
def mergeRateMonoR(f,M):
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
    return (3.7e6)*(f*0.85)**2*(M)**(-32/37)*((f*0.85)**2+0.005**2)**(-21/74)
#Minimum and maximum mass of primordial black hole 
Mmin=10;
Mmax=100;
#counter
l=0;
Nmax=150;# number of maximum size of the vector of  total merge rate.
G=np.zeros((Nmax,1)); # Power-law
K=np.zeros((Nmax,1)); #Log normal
D=np.zeros((Nmax,1)); #CriticalCollapse
fx=np.linspace(1e-4,1,Nmax) #Values of f that we have considered
#Loop to obtain the total merge rate for different values of f_PBH
for f in fx:
  def mergeRatePowerLawf(Mi,Mj):
       return mergeRatePLBlaw(Mi,Mj,f,I)
  def mergeRatelogNormalf(Mi,Mj):
        return mergeRatePLB(Mi,Mj,f,C)
  def mergeRateCollapsef(Mi,Mj):
        return mergeRatePLcc(Mi,Mj,f,H)      
  G[l]=din(mergeRatePowerLawf,Mmin,Mmax,Mmin,Mmax)[0]
  K[l]=din(mergeRatelogNormalf,Mmin,Mmax,Mmin,Mmax)[0]
  D[l]=din(mergeRateCollapsef,Mmin,Mmax,Mmin,Mmax)[0]
  l=l+1;
  print(l)
  
#Plot save
plt.plot(fx,mergeRateMonoR(fx,30),'b',label="Monochromatic") 
plt.plot(np.linspace(1e-4,1,Nmax),G,'k',label="Power-Law")
plt.plot(np.linspace(1e-4,1,Nmax),K,'C1',label="Lognormal")
plt.plot(np.linspace(1e-4,1,Nmax),D,'g',label="Critical Collapse")
plt.xlabel("$f_{PBH}=\Omega_{PBH}\: / \:\Omega_{DM}$")
plt.ylabel(" Merge Rate [$Gpc^{-3} yr^{-1}$]")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig("Plots/TotalFigure.eps")
plt.savefig("Plots/TotalFigure.pdf")

