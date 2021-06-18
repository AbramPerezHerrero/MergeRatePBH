# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 20:04:55 2021

@author: aph278
"""
from scipy import optimize
import MassFunction as MF
import numpy as np
import ReaderTxt as RT
from scipy.interpolate import interp1d

Nmax=30

v1=np.linspace(1e-4,1e-3,Nmax)
v2=np.linspace(1e-3,1e-2,Nmax)
v3=np.linspace(1e-2,1e-1,Nmax)
v4=np.linspace(1e-1,1,Nmax)
fx=np.concatenate((v1[0:29],v2[0:29],v3[0:29],v4))
#O2
Ta=0.428
McArti=np.loadtxt('19060800/19060800.txt',dtype=float)[:,0]
fArti=np.loadtxt('19060800/19060800.txt',dtype=float)[:,1]
R90=np.zeros((len(fArti)))
for i in np.arange(0,len(fArti),1):
 R90[i]=fArti[i]**(-3)*(-np.log(0.1)*3/(4*np.pi*Ta))

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

# def Monochromatic(M,I,M01):
#     '''
#     Mass function used which takes into account normalisation, Mc and 
#     alpha value. 

#     Parameters
#     ----------
#     M : FLOAT
#         Primordial black hole mass in solar units
#     I : FLOAT
#         Normalisation of the mass function taking into account the parameters.
#     delta1 : FLOAT
#         delta of the mass function that has been chosen.
#     alpha1: FLOAT
#         alpha parameter of the PL mass function
#     Returns
#     -------
#     FLOAT
#         The probability of the mass M.

#     '''
#     return MF.pseudomonochromatic(M,B=I,M0=M01)
# # This function represent the equation 3.58 of the work
# def mergeRate(Mi,Mj,f,N,M01):
#     # 'zeta' correpond to the clustering parameter. If 'zeta=1' correspond
#     # to non-clustering
#     zeta=1 
#     '''
#     Merge rate obtained in the work.

#     Parameters
#     ----------
#     Mi : FLOAT
#         Correspond to M1 in the work, are the PBH's mass in solar units.
#     Mj : FLOAT
#         Correspond to M2 in the work, are the PBH's mass in solar units.
#     f : FLOAT
#         The abundance of primordial black holes. It must be in the interval
#         [0,1]
#     N : FLOAT
#         Normalization of the mass function.
#     M01 : FLOAT
#         M0 parameter of the mass function that has been chosen.

#     Returns
#     -------
#     FLOAT
#         The merge rate of PBH with extended masss function. 

#     '''
#     # In this case the function used is monochromatic
#     function=Monochromatic
#     if Mj!=Mi:
#             return (3.7e6)*(f*0.85)**2*zeta**(16/37)*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*min(function(Mi,N,M01)/Mi,function(Mj,N,M01)/Mj)*(function(Mi,N,M01)/Mi+function(Mj,N,M01)/Mj)
#     else :
#             return (3.7e6)*(f*0.85)**2*zeta**(16/37)*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*(function(Mj,N,M01)/(2*Mj))*(function(Mj,N,M01)/(Mj))

Mmin=0
Mmax=100
a=0

fval=np.zeros((len(R90)))
A=np.zeros((len(fx))); 
for M01 in McArti:
    I=MF.normalization("monochromatic",M0=M01)
    l=0
    for f in fx:
        # def mergeRatef(Mj,Mi):
        #     return mergeRate(Mi,Mj,f,I,M01)   
        A[l]=mergeRateMono(f, M01) #dblquad(mergeRatef,Mmin,Mmax,Mmin,Mmax)[0]
        l=l+1
    inter=interp1d(fx, A, kind='linear')
    def solution(f):
        return inter(f)-R90[a]
    try:fval[a]=optimize.brentq(solution,1e-4,1) 
    except ValueError: fval[a]=1
    print(fval[a])
    a=a+1
namefile='MonochromaticHighMass'
data = np.column_stack([McArti,fval])
datafile_path = 'listfile/' + namefile+'.txt'
np.savetxt(datafile_path , data, fmt='%10.10f')

