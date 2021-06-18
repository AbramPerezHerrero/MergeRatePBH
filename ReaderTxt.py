# -*- coding: utf-8 -*-
"""
Program to open txt data.

@author: Abram Perez Herrero
"""
import numpy as np
def opentxt(file,separator,cutoff="yes"):
    '''
    Read the txt data in two colums.

    Parameters
    ----------
    file : STRING
        Name of the txt and his path.

    Returns
    -------
    list
        List with two data arrays.

    '''
    A=open(file,"r")
    lines=A.readlines()
    M=[]
    f=[]
    
    for x in lines:
           try:
               M.append(x.split(separator)[0])
               f.append(x.split(separator)[1])
           except IndexError:
               break
         
    A.close()
    m=0
    j=len(M)
    for i in np.arange(0,len(M),1): 
            M[i]=float(M[i])
            f[i]=float(f[i])
    for i in np.arange(0,int(len(M)/2),1):
                if f[i]==1 and cutoff=='yes':
                    m=i
    for i in np.arange(int(len(M)/2),len(M),1):
                if f[i]==1 and cutoff=='yes':
                    j=i  

    return [np.array(M)[m:j],np.array(f)[m:j]]