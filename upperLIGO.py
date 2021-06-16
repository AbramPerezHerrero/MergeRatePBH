# -*- coding: utf-8 -*-
"""
@author: Abram Perez Herrero
@based on: https://dcc.ligo.org/LIGO-P2000434/public.
@Reference:https://arxiv.org/pdf/2010.14533.pdf.
"""
import numpy as np
import deepdish as dd
def  upperLIGO(Data):
    '''
    Function that obtains the upper limit of the credibility interval of 
    the article https://arxiv.org/pdf/2010.14533.pdf. specifically for
    the Power-Law +Peak mass function.

    Returns
    -------
     Numpy array
        Array with the upper limit data provided by the LIGO/Virgo article.
        

    '''
    filenames = [
        Data+"/o1o2o3_mass_c_iid_mag_two_comp_iid_tilt_powerlaw_redshift_mass_data.h5",
        ]        
    peak_1 = 0
    _peak_1 = []    
    limits = [5,95] 
    ff=filenames[0]
    h = dd.io.load(ff)
    lines = h["lines"]
    print("Data analysed from {}".format(ff))
    _peak_1.append(max(np.percentile(lines["mass_1"], limits[1], axis=0)))
    peak_1 = max(peak_1, max(_peak_1))
    return np.percentile(lines["mass_1"], limits[1], axis=0)