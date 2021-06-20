'''
-This program serve to obtain the value of the maximum black hole abundance for 
 a Power-Law mass function and for values of \alpha=[4,2.7,2.2] 
 parameter.

@Return: A txt file with the values of  maximum f_{PBH} for each value of delta.
@author:Abram Pérez Herrero
@Date:14/06/21
'''


from time import time
from scipy.integrate import quad as din
import MassFunction as MF
import numpy as np
import upperLIGO as UL


# Variable 'v', represent a numpy array with all values of f_pbh which 
# we will consider.
Nmax=50; # Nmax, are the maximum number of point for each interval of  'v'.
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

def PL(M,I,delta1,alpha1):
    '''
    Mass function used which takes into account normalisation, Mc and 
    alpha value. 

    Parameters
    ----------
    M : FLOAT
        Primordial black hole mass in solar units
    I : FLOAT
        Normalisation of the mass function taking into account the parameters.
    delta1 : FLOAT
        delta of the mass function that has been chosen.
    alpha1: FLOAT
        alpha parameter of the PL mass function
    Returns
    -------
    FLOAT
        The probability of the mass M.

    '''
    return MF.powerLaw(M,A=I,alpha=alpha1,delta=delta1)

# This function represent the equation 3.58 of the work
def mergeRate(Mi,Mj,f,N,delta1,alpha1):
    # 'zeta' correpond to the clustering parameter. If 'zeta=1' correspond
    # to non-clustering
    zeta=1 
    '''
    Merge rate obtained in the work.

    Parameters
    ----------
    Mi : FLOAT
        Correspond to M1 in the work, are the PBH's mass in solar units.
    Mj : FLOAT
        Correspond to M2 in the work, are the PBH's mass in solar units.
    f : FLOAT
        The abundance of primordial black holes. It must be in the interval
        [0,1]
    N : FLOAT
        Normalization of the mass function.
    Mc1 : FLOAT
        Mc parameter of the mass function that has been chosen.
    beta1:
        beta parameter of the mass function
    Returns
    -------
    FLOAT
        The merge rate of PBH with extended masss function. 

    '''
    # In this case the function used is PL
    function=PL
    if Mj!=Mi:
            return (3.7e6)*(f*0.85)**2*zeta**(16/37)*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*min(function(Mi,N,delta1,alpha1)/Mi,function(Mj,N,delta1,alpha1)/Mj)*(function(Mi,N,delta1,alpha1)/Mi+function(Mj,N,delta1,alpha1)/Mj)
    else :
            return (3.7e6)*(f*0.85)**2*zeta**(16/37)*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*(function(Mj,N,delta1,alpha1)/(2*Mj))*(function(Mj,N,delta1,alpha1)/(Mj))



# This part corresponds to the ligo data, in order to obtain the upper limit 
# of dR/dm1. If an error occurs, check that the data is in the correct location.
# Only the data from the POWERL-LAW + PEAK graph has been considered.
# The name filedata correspond to the file in where the LV data are.
mass_1 = np.linspace(2, 100, 1000)
filedata='DataLVO3'
# 'upper' correspond to the upper limit of the credible interval of LV analysis
upper=UL.upperLIGO(filedata)
p=0 #The variable p is only a counter to show in terminal in which interation
    #are.
Mmin=3; #Mmin and Mmax are the cut_off in the integration of the merger rate.
Mmax=100;
# Array with the values of Mc
Mcvalues=np.concatenate((np.linspace(2,100,35)[0:34]
                         ,np.linspace(100,700,30)))
#The array fval is used to save the value of f_{pbh} for each value of Mc.
fval=np.zeros([len(Mcvalues)])

def check(vec,delta2,I,alpha1):
            '''
            Function that returns whether the function exceeds the limit in the 
            interval or group proposed by the vector or not.

            Parameters
            ----------
            v: Array
                Array with the interval of f in which this function is going
                to check.
            delta2 : FLOAT
                delta value of the mass function parameter.
            I : FLOAT
                Normalization of the mass function.
            alpha1 : FLOAT
                alpha parameter of the mass function.

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
                    return mergeRate(Mi,Mj,f,I,delta2,alpha1)   
                A[l]=din(mergeRatePLBf,Mmin,Mmax)[0]
                l=l+1;
            if(np.max(A[39:850]-upper[39:850])<=0):
                return False
            else:
                return True
#Names of the txt files where the data is saved            
names=['P-L4','P-L27','P-L22']
# A counter variable
h=0
# First loop for each value of parameter beta
for alpha1 in [4,2.7,2.2]:
    namefile=names[h]
    h=h+1
    # Second loop for each value of parameter Mc
    for delta1 in Mcvalues:
        print('$\\delta$:',delta1) # Only if you want to show in which step are analised
        # Obtain the normalization of the mass function
        I=MF.normalization("powerLaw",alpha=alpha1,delta=delta1)
        #time counter
        start_time = time()
        #step counter
        k=0
        j=30
        # It is check in which interval are the f_PBH maximum
        for i in np.arange(0,15,1):
            if check(v[i],delta1,I,alpha1):
                j=i
                print('f_{pbh} are between', v[j][0], 'and', v[j][Nmax-1]) 
                break;
        if(j==30):
            fval[p]=1
        else:  
            for f in v[j]:  #Third loop for each value of parameter f in the 
                            # interval
                # Array with the values of dR/dM1
                A=np.zeros((1000)); 
                # Counter
                l=0
                for Mi in mass_1:# Four loop for each value of M1 
                        def mergeRatef(Mj):
                            return mergeRate(Mi,Mj,f,I,delta1,alpha1)   
                        A[l]=din(mergeRatef,Mmin,Mmax)[0]
                        l=l+1;
                # It is compare with the upper limit in the range [6,95] Msun       
                if(np.max(A[39:850]-upper[39:850])<=0):       
                    fval[p]=f                
                else:
                    break; 
                k=k+1
                print('nº steps in f:',k)  # Only if you want to show 
                                           #in which step are analised  

                          
        print('nº step in Mc:',p)
        p=p+1            
        elapsed_time = time() - start_time
        print('time(s):',elapsed_time)
    # #Save the data in a txt document.   
    # data = np.column_stack([Mcvalues, fval])
    # datafile_path = 'listfile/' + namefile+'.txt'
    # np.savetxt(datafile_path , data, fmt='%10.10f')
    # fval=np.zeros([len(Mcvalues)])
    # p=0

    
