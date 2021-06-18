# -*- coding: utf-8 -*-
"""
Code for plot  the constraints of PBH assuming Monochromatic estimed.
@author: Abram Pérez Herrero
@Reference: https://github.com/bradkav/PBHbounds/
@Date:14/06/21
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import ReaderTxt as RT
from scipy.interpolate import interp1d
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
mpl.rcParams['legend.edgecolor'] = 'inherit'
plt.axhspan(1, 1.5, facecolor='grey', alpha=0.5)  
plt.ylim(1e-5, 1.5)
xmin = 1e-2
xmax = 1e3
plt.xlim(xmin, xmax)
plt.figure(figsize=(6,6))

# Load data, it is import that the data are dowload
Mclow=RT.opentxt("listfile/MonochromaticLowMass.txt"," ")[0]
flow=RT.opentxt("listfile/MonochromaticLowMass.txt"," ")[1]
Mc=RT.opentxt("listfile/Monochromatic.txt"," ")[0]
f=RT.opentxt("listfile/Monochromatic.txt"," ")[1]
fO3=RT.opentxt("listfile/MonochromaticLowMassO3.txt"," ")[1]
McO3=RT.opentxt("listfile/MonochromaticLowMassO3.txt"," ")[0]
fHigh=RT.opentxt("listfile/MonochromaticHighMass.txt"," ")[1]
MHigh=RT.opentxt("listfile/MonochromaticHighMass.txt"," ")[0]

finter=np.concatenate((flow[82:83],f[0:1]))
Mcinter=np.concatenate((Mclow[82:83],Mc[0:1]))

finter2=np.concatenate((fO3[len(fO3)-1:len(fO3)],f[0:1]))
Mcinter2=np.concatenate((McO3[len(fO3)-1:len(fO3)],Mc[0:1]))

finter3=np.concatenate((f[len(f)-1:len(f)],fHigh[0:1]))
Mcinter3=np.concatenate((Mc[len(f)-1:len(f)],MHigh[0:1]))
# finter=interp1d(Mcinter,finter,kind='linear')
# Mcinter=np.linspace(Mclow[82],Mc[0],10)



plt.fill_between(Mcinter2,finter2,1, color='red',alpha=0.2)
plt.plot(Mcinter2,finter2,'r--',label=r'Interpolation')

plt.fill_between(Mcinter3,finter3,1, color='red',alpha=0.2)
plt.plot(Mcinter3,finter3,'r--')

plt.fill_between(McO3,fO3,1, color='red',alpha=0.2)
plt.plot(McO3,fO3,'red')


plt.fill_between(Mc,f,1, color='red',alpha=0.2)
plt.plot(Mc,f,'r')

plt.fill_between(MHigh,fHigh,1, color='red',alpha=0.2)
plt.plot(MHigh,fHigh,'red',label=r'Result')

plt.legend(loc="lower left")
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.tick_bottom()
ax.xaxis.set_tick_params(pad=5)
ticks_minor = np.geomspace(1e-18, 1e4, 23)
ticks_minor = ticks_minor[(xmin <= ticks_minor) & (ticks_minor <= xmax)]
#print(ticks_minor)
ax.set_xticks(ticks_minor,minor=True)
ax.set_xticklabels([], minor=True)
plt.ylim(1e-5, 1.5)  
ax.set_xlabel(r'$M_\mathrm{0}$ [$M_\odot$]')
plt.ylabel(r'$f_\mathrm{PBH} = \Omega_\mathrm{PBH}/\Omega_\mathrm{DM}$')

ax_top = ax.twiny()
ax_top.xaxis.tick_top()
ax_top.set_xscale('log')
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xlabel(r'$M_\mathrm{PBH}$ [g]', labelpad=7)

g_to_Msun = 1/1.989e+33

g_ticks_minor = np.geomspace(1e15, 1e37, 23)
g_ticks_minor = g_ticks_minor[(xmin <= g_to_Msun*g_ticks_minor) & (g_to_Msun*g_ticks_minor <= xmax)]
g_ticks = g_ticks_minor


g_tick_labels = [r"$10^{" + str(int(np.log10(x))) +"}$" for x in g_ticks]

ax_top.set_xticks(g_ticks*g_to_Msun)
ax_top.set_xticklabels(g_tick_labels)
ax_top.xaxis.set_tick_params(pad=0)

ax_top.set_xticks(g_ticks_minor*g_to_Msun,minor=True)
ax_top.set_xticklabels([],minor=True)

plt.savefig("Plots/constraintsMonochromaticTotal.pdf", bbox_inches='tight')
    