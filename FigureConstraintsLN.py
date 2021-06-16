# -*- coding: utf-8 -*-
"""
Code for plot  the constraints of PBH assuming Log-Normal.
@author: Abram Pérez Herrero
@Reference: https://github.com/bradkav/PBHbounds/
@Date:14/06/21
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import ReaderTxt as RT
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
Mc1=RT.opentxt("listfile/LN1.txt",' ')[0]
f1=RT.opentxt("listfile/LN1.txt",' ')[1]

Mc06=RT.opentxt("listfile/LN06.txt",' ')[0]
f06=RT.opentxt("listfile/LN06.txt",' ')[1]

Mc04=RT.opentxt("listfile/LN04.txt",' ')[0]
f04=RT.opentxt("listfile/LN04.txt",' ')[1]



plt.fill_between(Mc1,f1,1, color='red',alpha=0.2)
plt.plot(Mc1,f1,'r',label="$\\sigma=1$")
plt.fill_between(Mc06,f06,1, color='darkblue',alpha=0.2)
plt.plot(Mc06,f06,'b',label="$\\sigma=0.6$")
plt.fill_between(Mc04,f04,1, color='mediumpurple',alpha=0.2)
plt.plot(Mc04,f04,'indigo',label="$\\sigma=0.4$")
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
ax.set_xlabel(r'$M_\mathrm{c}$ [$M_\odot$]')
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

plt.savefig("Plots/constrainsLN.pdf", bbox_inches='tight')
    