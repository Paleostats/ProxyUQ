#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 11:45:02 2022

@author: andreschristen

Program to analyze HP (Nicole) Cores
"""

from numpy import sum as npsum
from numpy import isnan

from matplotlib.pylab import plot, subplots, hist, close

from matplotlib.colors import TABLEAU_COLORS #CSS4_COLORS
colors = ["white"]+list(TABLEAU_COLORS)+["black"]

from ProxyUQ import ProxyUQ
from pyplotfuncmc import PlotSolnsMC


# Core ids.  HP1A has a problem
core_tag = ["1A", "1B", "1C", "1D", "2A", "2B", "2C", "2D", "3A", "3B", "3C", "3D"]

"""
for tag in core_tag:
    ProxyUQ( "HP/HP%s/HP%s" % (tag,tag), folder="Cores", p_list=None)

cont = {} #a container dictionary to hold all results
### Analyze all cores
for tag in core_tag:
    if tag == "1A":
        continue
    print("Core %s" % (tag), end="")
    rt = ProxyUQ( "HP/HP%s/HP%s" % (tag,tag), folder="Cores", p_list=[1], y_by=10)
    cont[tag] = rt
"""

### Plot all proxies in the same axis
fig, ax = subplots()
for i,tag in enumerate(core_tag):
    if tag == "1A":
        continue
    PlotSolnsMC( cont[tag]['dates'], cont[tag]['proxy']['Cdensity'], label=tag,\
                q=[0.1], fill=False, color=colors[i], med_col=colors[i], ax=ax)
ax.set_xlabel("y BP")
ax.set_ylabel("C density ($gr/cm^3$)")
ax.set_xlim((-80,700))
ax.legend(loc='upper right', ncol=3, fancybox=True, shadow=True)


### Plot the added C density histograms in a range of years 0 to 6 and 7 to 14
fig, ax = subplots(nrows=2,ncols=1,sharex=True)
for i,tag in enumerate(core_tag):
    if tag == "1A":
        continue
    print(tag)
    tmp = npsum(cont[tag]['proxy']['Cdensity'][:,:7], axis=1)
    ax[0].hist( tmp[~isnan(tmp)], color=colors[i], histtype='step', density=True)
    tmp = npsum(cont[tag]['proxy']['Cdensity'][:,7:14], axis=1)
    ax[1].hist( tmp[~isnan(tmp)], color=colors[i], histtype='step', density=True)
ax[1].set_xlabel("C density, acc. ($gr/cm^3$)")
ax[0].yaxis.set_label_position("right")
ax[0].set_ylabel("%.0f - %.0f, y BP" % (cont[tag]['dates'][0],cont[tag]['dates'][6]) )
ax[1].yaxis.set_label_position("right")
ax[1].set_ylabel("%.0f - %.0f, y BP" % (cont[tag]['dates'][7],cont[tag]['dates'][13]) )
ax[1].set_xlim((0.1,0.4))


### Plot the added C density histograms in a range of years 0 to 6 and 7 to 14
### In one single plot
fig, ax = subplots()
for i,tag in enumerate(core_tag):
    if tag == "1A":
        continue
    print(tag)
    tmp = npsum(cont[tag]['proxy']['Cdensity'][:,:7], axis=1)
    ax.hist( tmp[~isnan(tmp)], color=colors[i], histtype='step', density=True)
    tmp = npsum(cont[tag]['proxy']['Cdensity'][:,7:14], axis=1)
    ax.hist( tmp[~isnan(tmp)], color=colors[i], histtype='step', density=True, linestyle='dashed')
ax.set_xlabel("C density, acc. ($gr/cm^3$)")
ax.set_xlim((0.1,0.4))

### Now the histograms of differences
fig, ax = subplots(nrows=4,ncols=3)
for k,tag in enumerate(core_tag): #enumerate(core_tag):
    if tag == "1A":
        continue
    #print(tag)
    tmp1 = npsum(cont[tag]['proxy']['Cdensity'][:,:7], axis=1)
    tmp1 = tmp1[~isnan(tmp1)]
    tmp2 = npsum(cont[tag]['proxy']['Cdensity'][:,7:14], axis=1)
    tmp2 = tmp2[~isnan(tmp2)]
    if tmp1.size < tmp2.size:
        tmp = tmp1 - tmp2[:tmp1.size]
    else:
        tmp = tmp1[:tmp2.size] - tmp2
    i = k % 4
    j = k // 4
    ax[i,j].hist(tmp)
    ax[i,j].set_xlim((-0.2,0.2))
    ax[i,j].set_title(tag)
fig.tight_layout()







