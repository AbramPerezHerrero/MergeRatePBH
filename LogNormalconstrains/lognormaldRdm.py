
from time import time
from scipy.integrate import quad as din
import matplotlib.pyplot as plt
import MassFunction2 as MF2
import numpy as np
import deepdish as dd
'''
This program tries to obtain the value of the maximum black hole abundance for 
a lognormal mass function taking into account different parameters.
'''

# A function is created where the sigma value to be analysed is entered. 

def lognormal(M,I,mu):
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
    mu : FLOAT
    Mc of the mass function that has been chosen.

    Returns
    -------
    FLOAT
        The probability of the mass given.

    '''
    return MF2.logNormal(M,A=I,sigma=0.8,Mc=mu)
def mergeRatePLB(Mi,Mj,f,C,mu):
    '''
    Merge rate obtained in the work.

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
            return (3.7e6)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*min(function(Mi,C,mu)/Mi,function(Mj,C,mu)/Mj)*(function(Mi,C,mu)/Mi+function(Mj,C,mu)/Mj)
    else :
            return (3.7e6)*(f*0.85)**2*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*(function(Mj,C,mu)/(2*Mj))*(function(Mj,C,mu)/(Mj))

'''
This part corresponds to the ligo data, in order to obtain the upper limit 
of dR/dm1. If an error occurs, check that the data is in the correct location.
Only the data from the POWERL-LAW + PEAK graph has been considered.
'''
mass_1 = np.linspace(2, 100, 1000)
colours = [ "#DE8F05"]
filenames = [
    "o1o2o3_mass_c_iid_mag_two_comp_iid_tilt_powerlaw_redshift_mass_data.h5",
   ]
        
peak_1 = 0
_peak_1 = []

limits = [5,95]

mass_1 = np.linspace(2, 100, 1000)
mass_ratio = np.linspace(0.1, 1, 500)

plt.figure(figsize=(13,7))
n=0
ff=filenames[0]
h = dd.io.load(ff)
ppd = h["ppd"]
lines = h["lines"]
mass_1_ppd = np.trapz(ppd, mass_ratio, axis=0)
ass_ratio_ppd = np.trapz(ppd, mass_1, axis=-1)
print("Data analysed from {}".format(ff))
_peak_1.append(max(np.percentile(lines["mass_1"], limits[1], axis=0)))
peak_1 = max(peak_1, max(_peak_1))

'''
-The variable 'upper' correspond to the upper R90% of the figure 3 of 
https://arxiv.org/pdf/2010.14533.pdf.

-The array fval is used to save the value of f_{pbh} for each value of
 Mc.

(It will be analysed for a single sigma but for many MCs.)

-The variable p is only a counter to show in terminal in which interation
 are.

-nMC reprensent the number of Mc that we will analyse between Mcmax and Mclow 

-Mmin and Mmax are the cut_off in the integration of the merger rate.
'''
upper=np.percentile(lines["mass_1"], limits[1], axis=0)
fval=np.zeros([15,1])
p=0
nMc=15
Mcmax=75
Mclow=7
Mmin=3;
Mmax=100;
'''
This part works as follows:

-First, a loop is used for each value of Mc. Also, we  compute the 
corresponding normalisation of the mass function. 

-Then, another loop is started for each value of f, in an interval that is 
known that dR/dm obtained are greater than the upper limit given by LIGO. 
(as small as possible)

-The vector that is going to contain the value of dR/dm is created, 
which we call 'A'.

-A value l is created as a counter and so fill the vector A for each value
 of f.

- Finally, a last loop where each value of the mass M1 is considered. 
A new function of the merger rate is also defined where the value of the mass 
M1 has been inserted. Finally, it is integrated between the limits.

-After this is done, it is compared with the value of the upper limit taken 
 from LIGO/Virgo. Note a small detail in that the functions we are dealing with
 extend too much to small masses and this causes them to go out of function
 very quickly (as happens with the function proposed by them in orange in
 figure 3). Therefore, we remove the first 24 points and compare whether the
 difference is greater than zero, i.e. whether it has passed the cut-off point.
 The value of f is then stored.
 
'''
for mu in np.linspace(Mclow,Mcmax,nMc):
    I=MF2.normalization("LogNormal",sigma=0.8,Mc=mu)
    k=0
    for f in np.linspace(0.0005,0.0011,150):      
        A=np.zeros((1000,1));
        l=0
        for Mi in mass_1:
            def mergeRatePLBf(Mj):
                return mergeRatePLB(Mi,Mj,f,I,mu)   
            A[l]=din(mergeRatePLBf,Mmin,Mmax)[0]
            l=l+1;
        if(np.max(A[24:1000,0]-upper[24:1000])<=0):       
                fval[p]=f                
        else:
            break;
        
        k=k+1 
        print(k)        
    print(p)
    p=p+1

'''
When the program finish, for each value of Mc the maximum value of f has been
stored in fval.
'''
    