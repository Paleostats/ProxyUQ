#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 09:41:10 2022

@author: jac
"""

from numpy import loadtxt, arange, zeros, cumsum, where, ceil, floor, array
from scipy.interpolate import interp1d

from matplotlib.pylab import plot, subplots, hist
from pandas import read_csv

from pyplotfuncmc import PlotSolnsMC



def ProxyUQ( core, p_list=[], y_by=10):
    """Produce a posterio MC sample of a proxy, from the MC output of
       an age-depth model, e.g. Bacon.
       core - Core prefix, e.g. "AtmosphericRivers/LL14", from this the
       file "AtmosphericRivers/LL14_settings.txt" and "AtmosphericRivers/LL14_proxies.csv" 
       will be read and deduce the .out file to "AtmosphericRivers/LL14_60.out"
       
       p_list - of empty, attempt to read the fikes and provied basic info.
              -otherwise, give the list of proxies to analyse, e.g. p_list=[1,5] 

       y_by - step for the age grid.
       
       returns a list with the age grid and the MC sample for each proxy analyzed
       
    """ 
    
    ### Read the first three lines of settings
    d_min, d_max, d_by = loadtxt("%s_settings.txt" % (core), max_rows=3)
    K1 = int(floor((d_max-d_min)/d_by))+1
    out = loadtxt( "%s_%d.out" % (core,K1), delimiter=' ')
    T, tmp = out.shape
    K = tmp-3
    
    ### Read the proxie file and, for the moment, change missing values to 0
    proxy = read_csv("%s_proxies.csv" % (core))
    proxy = proxy.fillna(0)
    proxy_names = list(proxy.columns)
    proxy=array(proxy)
    if p_list == []:
        print("\n%s_%d.out" % (core,K1))
        print("d_min=%f, d_max=%f, d_by=%f, K=%d, T=%d" % (d_min,d_max,d_by,K,T))
        print("Proxies: ", proxy_names, "\n")
        return

    k_proxy = len(p_list)
    ### age grid for age-depth models
    depths = arange( d_min, d_max+d_by, step=d_by)
    
    ### Array to hold the age-depth models
    ### Sample size T, num of rows, and $K$ number of ages (Bacon sections)
    solns = zeros((T,K))
    
    for t in range(T): #make age depth models
        solns[t,:] = out[t,0] + cumsum( d_by*out[t,1:(K+1)])
    
    ### age grid for the proxies
    CalA = arange( ceil(solns.min()), floor(solns.max()), step=y_by)
    ### Array to hold the proxy MC samples
    ### Sample size, T rows, and CalA.size number of ages.
    solns_p = zeros((T,CalA.size))
    
    ### Plot the age-depth models with quantiles
    fig, ax = subplots( nrows=k_proxy+1, ncols=1, sharex=True) #, sharey=True)
    PlotSolnsMC( depths, solns, ax=ax[0]) #Default quantiles 0.1,0.25,0.5,0.75,0.9
    ax[-1].set_xlabel("Depth (cm)")
    ax[0].set_ylabel("Age (y BP)")
    
    fig_p, ax_p = subplots( nrows=k_proxy, ncols=1, sharex=True, squeeze=False)
    ax_p[-1,0].set_xlabel("Age (y BP)")
    
    rt = [CalA]
    ### Iterate through the proxy list
    for ax_n,p in enumerate(p_list):
        ### Interpolate the proxy
        inter_p = interp1d( proxy[:,0], proxy[:,p]) 
        print("\nProxy (%s)" % (proxy_names[p]))
        ### Itearete through the calendar ages
        for i,g in enumerate(CalA):
            ### For each year, iterate through all the age-depth models
            print(g, end=' ')
            for t in range(T):
                ### Here we find the inverse d
                tmp = where(solns[t,:] < g)[0]
                if len(tmp) == 0:
                    indx = K-2
                else:
                    indx = tmp[-1]
                d = depths[indx] + (g-solns[t,indx])/out[t,indx+2]
                ### Find the corresponding proxy value
                if (proxy[0,0] <= d) and (d <= proxy[-1,0] ):
                    solns_p[t,i] = inter_p(d) #add error here for proxies with error!
                else:
                    solns_p[t,i] = None
    
        ax[ax_n+1].plot( proxy[:,0], proxy[:,p], '-')
        ax[ax_n+1].set_ylabel("Proxy (%s)" % (proxy_names[p]))
        
        ### Plot the proxy ate cal age with quentiles
        PlotSolnsMC( CalA, solns_p, ax=ax_p[ax_n,0])
        ax_p[ax_n,0].set_ylabel("Proxy (%s)" % (proxy_names[p]))
        rt += [solns_p.copy()]
    
    fig.tight_layout()
    fig_p.tight_layout()
    return rt


if __name__ == "__main__":
    
    ### Read the files and provide general info 
    ProxyUQ( "HP1C/HP1C", p_list=[])
    ProxyUQ( "Auassat/Auassat", p_list=[])
    ProxyUQ( "AtmosphericRivers/LL14", p_list=[])

    """
    ### Example:
    ### Analyze proxies 2 and 3.  Note, start at 1, 0 is the depth column.
    CalA, solns_AC, solns_CAR = ProxyUQ( "Auassat/Auassat", p_list=[2,3], y_by=10)
    ### You can make the nice plot with:
    PlotSolnsMC( CalA, solns_AC)
    ### Or make a histogram at a particula age:
    fig, ax = subplots()
    ax.hist(solns_AC[25,:]) #Proxy at age CalA[25]
    ax.set_xlabel("Acc. Carbon")
    ax.set_title("%.0f, y BP" % (CalA[25]))
    """
    
    

