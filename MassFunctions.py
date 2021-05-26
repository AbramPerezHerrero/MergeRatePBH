'''
CREATED BY ABRAM PÉREZ HERRERO 
VERSION 1 , 26/05/21
THIS PACKAGE INCLUDE THE PBH MASS FUNCTIONS:
    -powerLaw
    -logNormal
    -criticalCollapse
    -powerLawBroken
    -powerLawPeak
    -multipeak
    -Monocrhomatic
'''

import numpy as np
from scipy.integrate import quad as inte

def funf(delta,M,Mmin):
    return np.exp((delta/(M-Mmin))+(delta/(M-Mmin-delta)))
def S(M,Mmin,delta):
  if M<Mmin:
    return 0.0
  elif Mmin<=M and M<Mmin+delta:
    return (funf(delta,M,Mmin)+1)**(-1)
  elif M>(Mmin+delta):
    return 1
def Gaussian(M,sigma,mu):
  return (1/np.sqrt(2*np.pi*sigma**2))*np.exp(-(M-mu)**2/(2*sigma**2))
def Powerlaw2(M,Mmax,alpha):
 if M<Mmax:
     return M**(-alpha)
 else:
     return 0

def powerLaw(M,Mmin=4.59,alpha=1.6,delta=4.85):
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
    def powerLawnormal(M):
        #return((alpha-1)/(M))*(M/Mc)**(-alpha)*S(M,Mmin,delta)
        return M**(-alpha)*S(M,Mmin,delta)
    return powerLawnormal(M)/inte(powerLawnormal,0,100)[0]

def logNormal(M,sigma=0.6,Mc=15):
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
    def logNormalnormal(M):
      return (1/(np.sqrt(2*np.pi)*sigma*M))*np.exp(-(np.log10(M/Mc))**2/(2*sigma**2))
    return logNormalnormal(M)/inte(logNormalnormal,0,100)[0]

def criticalCollapse(M,beta=2.85,Mc=40):
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
    def criticalCollapsenormal(M):
     return M**beta*np.exp(-(M/Mc)**beta)
    return criticalCollapsenormal(M)/inte(criticalCollapsenormal,0,100)[0]
def powerLawBroken(M,Mmin=4.59,Mmax=86.22,b=0.43,delta=4.82,alpha1=1.58,alpha2=5.59):
    '''
    BROKEN POWER LAW MASS DISTRIBUTION FUNCTION. IF b=1 AND delta=0 THEN 
    IS EQUAL TO THE TRUNCATED POWER LAW MASS DISTRIBUTION FUNCTION.
    Parameters
    ----------
    M : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    Mmin : FLOAT
        MINIMUM MASS OF PBHS IN SOLAR MASS UNIT.IT IS 
        USED AS LOWER MASS CUT-OFF.
    Mmax : FLOAT
        MAXIMUM MASS OF PBHS IN SOLAR MASS UNIT. IT IS 
        USED AS HIGH MASS CUT-OFF.
    b : FLOAT
        FRACTION AT WHICH THE MASS DISTRIBUTION BREAKS. IF b=1 THE BREAK MASS
        ARE EQUAL TO THE MAXIMUM PBH MASS Mmax.
    delta : FLOAT
        TUNING RANGE FOR LOWER MASSES IN SOLAR MASS UNIT.
    alpha1 : FLOAT
        POWER-LAW INDEX USING IN THE RANGE OF MASS LOWER THAN THE BREAK MASS.
        
    alpha2 : FLOAT
        POWER-LAW INDEX USING IN THE RANGE OF MASS HIGHER THAN THE BREAK MASS.

    Returns
    -------
    FLOAT
        THE PROBABILITY OF SUCH MASS M TAKING INTO ACCOUNT THE BROKEN POWER LAW
        MASS DISTRIBUTION.

    '''
    Mbreak=Mmin+(b)*(Mmax-Mmin)
    def powerLawBrokennormal(M):
        if Mmin<Mmax and M<Mbreak:
            return M**(-alpha1)*S(M,Mmin,delta)
        elif Mbreak<M and M<Mmax:
            return M**(-alpha2)*S(M,Mmin,delta)
        else:
            return 0
    return powerLawBrokennormal(M)/inte(powerLawBrokennormal,0,100)[0]

def powerLawPeak(M,Mmin=4.59,Mmax=86.22,peak=0.1,alpha=2.63,sigma=5.69,mu=33.07,delta=4.82):
    '''
    POWER LAW+PEAK MASS DISTRIBUTION FUNCTION

    Parameters
    ----------
    M : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    Mmin : FLOAT
        MINIMUM MASS OF PBHS IN SOLAR MASS UNIT.IT IS 
        USED AS LOWER MASS CUT-OFF.
    Mmax : FLOAT
        MAXIMUM MASS OF PBHS IN SOLAR MASS UNIT. IT IS 
        USED AS HIGH MASS CUT-OFF.
    peak : FLOAT
        RATIO OF THE NUMBER OF BLACK HOLES AT THE GAUSSIAN PEAK 
        COMPARED TO THE REMAINDER OF THE FUNCTION. IT MUST BE BETWEEN 0 TO 1.
    alpha : FLOAT
        POWER LAW INDEX.
    sigma : FLOAT
        STANDAR DEVIATION OF THE GAUSSIAN PEAK.
    mu : FLOAT
        EXPECTED VALUE OF THE GAUSSIAN PEAK.
    delta : FLOAT
        TUNING RANGE FOR LOWER MASSES IN SOLAR MASS UNIT.

    Returns
    -------
    FLOAT
        THE PROBABILITY OF SUCH MASS M TAKING INTO ACCOUNT THE POWER-LAW +PEAK
        MASS DISTRIBUTION.

    '''
    def powerlaw2nomal(M):
        return (Powerlaw2(M,Mmax,alpha))*S(M,Mmin,delta)
    def peakfun(M):
     return Gaussian(M,sigma,mu)*S(M,Mmin,delta)
    return (1-peak)*(powerlaw2nomal(M)/inte(powerlaw2nomal,0,100)[0]+peak*peakfun(M)/inte(peakfun,0,100)[0])

def multipeak(M,Mmin=4.59,Mmax=86.22,alpha=2.63,peak1=0.09,peak2=0.92,sigma1=5.69,mu1=33.39,sigma2=6.76,mu2=68.03,delta=4.82):
    '''
    MULTI-PEAK MASS DISTRIBUTION FUNCTION

    Parameters
    ----------
    M : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    Mmin : FLOAT
        MINIMUM MASS OF PBHS IN SOLAR MASS UNIT.IT IS 
        USED AS LOWER MASS CUT-OFF.
    Mmax : FLOAT
        MAXIMUM MASS OF PBHS IN SOLAR MASS UNIT. IT IS 
        USED AS HIGH MASS CUT-OFF.
    alpha : FLOAT
        POWER LAW INDEX.
    peak1 : FLOAT
        RATIO OF THE NUMBER OF BLACK HOLES AT THE GAUSSIAN  FIRST PEAK 
        COMPARED TO THE REMAINDER OF THE FUNCTION. IT MUST BE BETWEEN 0 TO 1.
    peak2 : FLOAT
        RATIO OF THE NUMBER OF BLACK HOLES AT THE GAUSSIAN SECOND PEAK 
        COMPARED TO THE REMAINDER OF THE FUNCTION. IT MUST BE BETWEEN 0 TO 1.
    sigma1 : FLOAT
        STANDAR DEVIATION OF THE FIRST GAUSSIAN PEAK.
    mu1 : FLOAT
        EXPECTED VALUE OF THE  FIRST GAUSSIAN PEAK.
    sigma2 : FLOAT
        STANDAR DEVIATION OF THE SECOND GAUSSIAN PEAK.
    mu2 : FLOAT
        EXPECTED VALUE OF THE SECOND GAUSSIAN PEAK.
    delta : FLOAT
        TUNING RANGE FOR LOWER MASSES IN SOLAR MASS UNIT.
    Returns
    -------
    FLOAT
        THE PROBABILITY OF SUCH MASS M TAKING INTO ACCOUNT THE MULTI-PEAK
        MASS DISTRIBUTION.
    '''
    def powerlaw2normal(M):
        return Powerlaw2(M,Mmax,alpha)*S(M,Mmin,delta)
    def peak1fun(M):
        return Gaussian(M,sigma1,mu1)*S(M,Mmin,delta)
    def peak2fun(M):
        return Gaussian(M,sigma2,mu2)*S(M,Mmin,delta)
    
    return ((1.0-peak1)*powerlaw2normal(M)/inte(powerlaw2normal,0,100)[0]+peak1*peak2*peak1fun(M)/inte(peak1fun,0,100)[0]+peak1*(1.0-peak2)*peak2fun(M)/inte(peak2fun,0,100)[0])
def monocrhomatic(M,M0=30.0):
    '''
    MONOCHROMATIC MASS FUNCTION.

    Parameters
    ----------
    M : FLOAT
        PRIMORDIAL BLACK HOLE MASS IN SOLAR MASS UNIT.
    M0 : FLOAT
       SINGLE MASS OF PRIMORDIAL BLACK HOLES IN SOLAR MASS UNIT.

    Returns
    -------
    FLOAT
        THE PROBABILITY OF SUCH MASS M TAKING INTO ACCOUNT THE MONOCHROMATIC
        MASS FUNCTION.

    '''
    if(M<M0+0.1 and M>M0-0.1):
        return 1;
    else:
        return 0;

    