'''
THIS PACKAGE INCLUDE THE PBH MASS FUNCTIONS AND HIS NORMALIZATION:
    -normalization
    -powerLaw
    -logNormal
    -criticalCollapse

@author: Abram Pérez Herero
@Date: 14/06/2021

'''
import numpy as np
from scipy.integrate import quad as inte
# This functions are described in the work:
#    -g is the exponential function to obtain delta
#    -Gaussian is the typical gaussian distribution function
# They are use to describe the others mass function.
def funf(delta,M,Mmin):
    return np.exp((delta/(M-Mmin))+(delta/(M-Mmin-delta)))
def S(M,Mmin,delta):
  if M<Mmin:
    return 0.0
  elif Mmin<M and M<(Mmin+delta):
    return (funf(delta,M,Mmin)+1)**(-1)
  elif M>(Mmin+delta):
    return 1
def Gaussian(M,sigma,mu):
  return (1/np.sqrt(2*np.pi*sigma**2))*np.exp(-(M-mu)**2/(2*sigma**2))


def normalization(x,Mmin=3,Mmax=np.inf,Mc=15,alpha=2.63,delta=4.82,beta=2.85,
                  sigma=0.6,mu3=20,M0=30):
    '''
    Function that provides the normalisation value for each of the functions.
    Parameters
    ----------
    x : STRING
        STRING WITH NAME OF THE FUNCTION THAT YOU WANT TO NORMALIZATE. IT 
        COULD BE 'powerLaw', 'criticalCollapse', 'LogNormal or 'monochromatic'.
    Mmin : FLOAT
        MINIMUM MASS OF PBHS IN SOLAR MASS UNIT.IT IS 
        USED AS LOWER MASS CUT-OFF.
    Mmax : FLOAT
        MAXIMUM MASS OF PBHS IN SOLAR MASS UNIT.IT IS 
        USED AS LOWER MASS CUT-OFF.
    alpha : FLOAT
        POWER LAW SPETRAL INDEX.
    delta : FLOAT
        TUNING RANGE FOR LOWER MASSES IN SOLAR MASS UNIT.
    sigma : TYPE, optional
        DESCRIPTION. The default is 0.6.
    Mc : FLOAT
        MASS SCALE FACTOR IN SOLAR UNIT.
    beta : FLOAT
       EXPONENT THAT GIVES FORM TO THE FUNCTION.THE SMALLER IT IS THE WIDER 
       THE FUNCTION BECOMES .

    Returns
    -------
    The normalization of the funciton mass which the parameters given
    '''
     
    if x=="powerLaw":
        def powerLawnormal(M):
                 if( M>=Mmin):
                     return M**(-alpha)*S(M,Mmin,delta)
                 else:
                     return 0
        return inte(powerLawnormal,Mmin,np.inf)[0]
    elif  x=="LogNormal":
        def logNormalnormal(M):
                return (1/(np.sqrt(2*np.pi)*sigma*M))*np.exp(-(np.log(M/Mc))**2/(2*sigma**2))
        return inte(logNormalnormal,0,np.inf)[0]
    elif x=="criticalCollapse":
        def criticalCollapsenormal(M):
                if( M>Mmin):
                    return M**beta*np.exp(-(M/Mc)**beta)
                else:                    
                    return 0
        return inte(criticalCollapsenormal,0,np.inf)[0]
    elif  x=="monochromatic":
        def logNormalnormal(M):
                return (1/(np.sqrt(2*np.pi)*sigma*M))*np.exp(-(np.log(M/Mc))**2/(2*sigma**2))
        return inte(logNormalnormal,0,100)[0]

def powerLaw(M,Mmin=3,alpha=1.6,delta=4.85,A=1,Mmax=90):
    '''
    POWER-LAW MASS DISTRIBUTION FUNCTION.    
    Parameters
    ----------
    M : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    Mmin : FLOAT
        MINIMUM MASS OF PBHS IN SOLAR MASS UNIT.IT IS 
        USED AS LOWER MASS CUT-OFF.
    alpha : FLOAT
        POWER LAW SPETRAL INDEX.
    Mc : FLOAT
        MASS SCALE FACTOR IN SOLAR UNIT.
    delta : FLOAT
        TUNING RANGE FOR LOWER MASSES IN SOLAR MASS UNIT.

    Returns
    -------
    FLOAT
        THE PROBABILITY OF SUCH MASS M TAKING INTO ACCOUNT THE POWER-LAW
        MASS DISTRIBUTION.

    '''
    if(M>=Mmin):
      return M**(-alpha)*S(M,Mmin,delta)/float(A)
    else:
      return 0

def logNormal(M,sigma=0.6,Mc=15,A=1,Mmin=3,Mmax=90):
    '''
    LOG-NORMAL MASS DISTRIBUTION FUNCTION.
    Parameters
    ----------
    M : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    sigma : FLOAT
        THE STANDAR DEVIATION OF THE MASS DISTRIBUTION.
    Mc : FLOAT
        THE MEAN OF THE MASS DISTRIBUTION IN SOLAR UNIT.

    Returns
    -------
    FLOAT
        THE PROBABILITY OF SUCH MASS M TAKING INTO ACCOUNT THE LOG-NORMAL
        MASS DISTRIBUTION.

    '''

    return (1/(np.sqrt(2*np.pi)*sigma*M))*np.exp(-(np.log(M/Mc))**2/(2*sigma**2))/A


def criticalCollapse(M,Mmax=90,Mmin=3,beta=2.85,Mc=40,A=1):
    '''
    CRITICAL COLLAPSE MASS DISTRIBUTION FUNCTION.
    Parameters
    ----------
    M : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    beta : FLOAT
       EXPONENT THAT GIVES FORM TO THE FUNCTION.THE SMALLER IT IS THE WIDER 
       THE FUNCTION BECOMES .
    Mc : FLOAT
        EXPONENTIAL HIGH MASS CUT-OFF IN SOLAR MASS UNIT.

    Returns
    -------
    FLOAT
        THE PROBABILITY OF SUCH MASS M TAKING INTO ACCOUNT THE CRITICAL COLLAPSE
        MASS DISTRIBUTION.

    '''
    if( M>Mmin):
        return M**beta*np.exp(-(M/Mc)**beta)/A
    else:
       return 0

def pseudomonochromatic(M,B,M0):
    '''
    PSEUDO-MONOCHROMATIC MASS FUNCTION. IT IS USED THE GAUSSIAN FUNCTION
    WITH A LOW SIGMA TO REPRESENT A DIRAC DELTA.

    Parameters
    ----------
    M : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    M0 : FLOAT
        UNIC MASS OF PRIMORDIAL BLACK HOLE IN SOLAR MASS UNIT.

    Returns
    -------
    FLOAT
        THE PROBABILITY OF SUCH MASS M TAKING INTO ACCOUNT THE 
        MONOCHROMATIC MASS DISTRIBUTION.

    '''
    sigma=0.05
    return (1/(np.sqrt(2*np.pi)*sigma*M))*np.exp(-(np.log(M/M0))**2/(2*sigma**2))/B