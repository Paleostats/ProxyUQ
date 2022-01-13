#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 11:00:26 2022

@author: jac
"""

from numpy import array

from pandas import read_csv, unique

from ProxyUQ import ProxyUQ
from pyplotfuncmc import PlotSolnsMC


cat = read_csv("AtmosphericRivers/AR_Cat.csv")

ARivers = []
for y in unique(cat['Year']):
    tmp = cat.loc[cat['Year'] == y]
    intensity = 0
    for index, row in tmp.iterrows():
        #intensity += row['Hours']*row['Category']
        intensity += row['IWV max']
    #print("%d - %f" % (y,intensity))
    ARivers += [[y,intensity]]
ARivers = array(ARivers)




#rt = ProxyUQ( "AtmosphericRivers/LL14", p_list=[1,2,3,4,5,6], y_by=1)
ax, quan = PlotSolnsMC( rt[0], rt[4])
ax.plot( 1950-ARivers[:,0], ARivers[:,1]/5000+0.1, '-')
ax.set_ylabel("Al/Si")
ax.set_xlabel("y BP")

