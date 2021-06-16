from time import time
from scipy.integrate import quad as din
import matplotlib.pyplot as plt
import MassFunction4 as MF
import numpy as np
import deepdish as dd
'''
This program tries to obtain the value of the maximum black hole abundance for 
a lognormal mass function taking into account different parameters.

v represent a vector with all values of f_pbh which we will consider.
They are separated by orders of magnitude.
'''
Nmax=50;
v5=np.linspace(1e-5,3e-5,Nmax)
v55=np.linspace(3e-5,6e-5,Nmax)
v555=np.linspace(6e-5,1e-4,Nmax)
v1=np.linspace(1e-4,3e-4,Nmax)
v11=np.linspace(3e-4,6e-4,Nmax)
v111=np.linspace(6e-4,1e-3,Nmax)
v2=np.linspace(1e-3,3e-3,Nmax)
v22=np.linspace(3e-3,6e-3,Nmax)
v222=np.linspace(6e-3,1e-2,Nmax)
v3=np.linspace(1e-2,3e-2,Nmax)
v33=np.linspace(3e-2,6e-2,Nmax)
v333=np.linspace(6e-2,1e-1,Nmax)
v4=np.linspace(1e-1,3e-1,Nmax)
v44=np.linspace(3e-1,6e-1,Nmax)
v444=np.linspace(6e-1,1,Nmax)
v=[v5,v55,v555,v1,v11,v111,v2,v22,v222,v3,v33,v333,v4,v44,v444]
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
    return MF.logNormal(M,A=I,sigma=0.6,Mc=mu)
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

-Mcvalues, are the array of values of the paramet Mc to analyse.
'''
upper=np.percentile(lines["mass_1"], limits[1], axis=0)
p=0
Mcvalues=np.concatenate((np.linspace(0.75,2,20)[0:20-1],np.linspace(2,450,50)[0:49],np.linspace(450,700,20)))
fval=np.zeros([len(Mcvalues)])
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
 figure 3). Therefore, we remove the first 82 points and compare whether the
 difference is greater than zero, i.e. whether it has passed the cut-off point.
 The value of f is then stored.
 
'''

def check(vec,mu2,I):
            '''
            Function that returns whether the function exceeds the limit in the 
            interval or group proposed by the vector or not.  For example,if the
            merge rate touch the upper limit in 2e-1 only return true in the 
            last group [1e-1,1] (See vector v)

            Parameters
            ----------
            v: Array
                Array with the interval of f in which this function is going
                to check.

             Returns
             -------
            Boolean
              Return True if the intersection are in this interval or False if 
              not.

             '''
            f=vec[Nmax-1]
            A=np.zeros((1000));
            l=0        
            for Mi in mass_1:
                def mergeRatePLBf(Mj):
                    return mergeRatePLB(Mi,Mj,f,I,mu2)   
                A[l]=din(mergeRatePLBf,3,100)[0]
                l=l+1;
            if(np.max(A[39:850]-upper[39:850])<=0):
                return False
            else:
                return True

                        
for mu in Mcvalues:
    print('Mc:',mu)
    I=MF.normalization("LogNormal",sigma=0.6,Mc=mu)
    k=0
    start_time = time()    
    for i in np.arange(0,16,1):
        if check(v[i],mu,I):
            j=i
            print('f_{pbh} are between', v[j][0], 'and', v[j][Nmax-1]) 
            break;
    for f in v[j]:      
        A=np.zeros((1000));
        l=0
        for Mi in mass_1:
            def mergeRatePLBf(Mj):
                return mergeRatePLB(Mi,Mj,f,I,mu)   
            A[l]=din(mergeRatePLBf,3,100)[0]
            l=l+1;
        if(np.max(A[39:850]-upper[39:850])<=0):       
                fval[p]=f                
        else:
            break;
        k=k+1 
        print("nº step of f:",k)        
    print("nº step of mu",p)
    p=p+1            
    elapsed_time = time() - start_time
    print(elapsed_time)


'''
When the program finish, for each value of Mc the maximum value of f has been
stored in fval.
'''
    