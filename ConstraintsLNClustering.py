'''
-This program serve to obtain the value of the maximum black hole abundance for 
 a Log-Normal mass function and for value of \sigma=[0.6] 
 parameter and parameter clustering [10,100]

@Return: A txt file with the values of  maximum f_{PBH} for each value of Mc.
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

def logNormal(M,I,Mc1,sigma1):
    '''
    Mass function used which takes into account normalisation, Mc and 
    beta value. 

    Parameters
    ----------
    M : FLOAT
        Primordial black hole mass in solar units
    I : FLOAT
        Normalisation of the mass function taking into account the parameters.
    Mc1 : FLOAT
        Mc of the mass function that has been chosen.
    sigma1: FLOAT
        sigma parameter of the LN mass function
    Returns
    -------
    FLOAT
        The probability of the mass M.

    '''
    return MF.logNormal(M,A=I,sigma=sigma1,Mc=Mc1)

# This function represent the equation 3.58 of the work
def mergeRate(Mi,Mj,f,N,Mc1,sigma1,zeta):
    # 'zeta' correpond to the clustering parameter. If 'zeta=1' correspond
    # to non-clustering
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
    sigma1:
        sigma parameter of the mass function
    zeta1:
       clustering parameter of the mass function
    Returns
    -------
    FLOAT
        The merge rate of PBH with extended masss function. 

    '''
    # In this case the function used is lognormal
    function=logNormal
    if Mj!=Mi:
            return (3.7e6)*(f*0.85)**2*zeta**(16/37)*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*min(function(Mi,N,Mc1,sigma1)/Mi,function(Mj,N,Mc1,sigma1)/Mj)*(function(Mi,N,Mc1,sigma1)/Mi+function(Mj,N,Mc1,sigma1)/Mj)
    else :
            return (3.7e6)*(f*0.85)**2*zeta**(16/37)*((f*0.85)**2+(0.005)**2)**(-21/74)*(Mi+Mj)**(36/37)*(Mi*Mj)**(3/37)*(function(Mj,N,Mc1,sigma1)/(2*Mj))*(function(Mj,N,Mc1,sigma1)/(Mj))



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
# #Mcvalues=np.concatenate((np.linspace(0.1,1,25)[0:25-1],
#                          np.linspace(1,95,35)[0:34],np.linspace(95,500,25)[0:24],
#                          np.linspace(450,700,20)))
    
Mcvalues=np.concatenate((np.linspace(0.1,7,20)[0:19],
                         np.linspace(7,95,16)[0:15],np.linspace(95,500,15)[0:14],
                         np.linspace(450,900,20)))    
#The array fval is used to save the value of f_{pbh} for each value of Mc.
fval=np.zeros([len(Mcvalues)])

def check(vec,Mc2,I,sigma1,zeta1):
            '''
            Function that returns whether the function exceeds the limit in the 
            interval or group proposed by the vector or not.

            Parameters
            ----------
            v: Array
                Array with the interval of f in which this function is going
                to check.
            Mc2 : FLOAT
                Mc value of the mass function parameter.
            I : FLOAT
                Normalization of the mass function.
            sigma1 : FLOAT
                Sigma parameter of the mass function.
            zeta1:
                clustering parameter of the mass function
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
                    return mergeRate(Mi,Mj,f,I,Mc2,sigma1,zeta1)   
                A[l]=din(mergeRatePLBf,Mmin,Mmax)[0]
                l=l+1;
            if(np.max(A[39:850]-upper[39:850])<=0):
                return False
            else:
                return True
#Names of the txt files where the data is saved            
names=['LN06Cl10','LN06Cl100']
# A counter variable
h=0
# First loop for each value of parameter sigma
for zeta1 in [10,100]:
    namefile=names[h]
    h=h+1
    # Second loop for each value of parameter Mc
    for Mc1 in Mcvalues:
        print('Mc:',Mc1) # Only if you want to show in which step are analised
        # Obtain the normalization of the mass function
        I=MF.normalization("LogNormal",sigma=0.6,Mc=Mc1)
        #time counter
        start_time = time()
        #step counter
        k=0
        j=30
        # It is check in which interval are the f_PBH maximum
        for i in np.arange(0,15,1):
            if check(v[i],Mc1,I,0.6,zeta1):
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
                            return mergeRate(Mi,Mj,f,I,Mc1,0.6,zeta1)   
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

                          
        print('nº step in mu:',p)
        p=p+1            
        elapsed_time = time() - start_time
        print('time(s):',elapsed_time)
    #Save the data in a txt document.   
    data = np.column_stack([Mcvalues, fval])
    datafile_path = 'listfile/' + namefile+'.txt'
    np.savetxt(datafile_path , data, fmt='%10.10f')
    fval=np.zeros([len(Mcvalues)])
    p=0

    
