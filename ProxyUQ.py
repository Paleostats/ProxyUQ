#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 09:41:10 2022

@author: jac
"""

from numpy import loadtxt, arange, zeros, cumsum, where, ceil, floor
from scipy.interpolate import interp1d

from matplotlib.pylab import plot, subplots, hist

from pyplotfuncmc import PlotSolnsMC



def ProxyUQ( core, p_list=[], y_by=10):
    
    d_min, d_max, d_by = loadtxt("%s_settings.txt" % (core), max_rows=3)
    tmp = int((d_max-d_min)/d_by)+1
    out = loadtxt( "%s_%d.out" % (core,tmp), delimiter=' ')
    T, tmp = out.shape
    K = tmp-3
    
    proxy = loadtxt("%s_proxies.csv" % (core), delimiter=',', skiprows=1)
    f = open("%s_proxies.csv" % (core))
    header = f.readline()[:-1] #Read the header without the /n
    f.close()
    proxy_names = header.split(',')
    if p_list == []:
        print("\n%s_%d.out" % (core,tmp))
        print("d_min=%f, d_max=%f, d_by=%f, K=%d, T=%d" % (d_min,d_max,d_by,K,T))
        print("Proxies: ", proxy_names, "\n")
        return

    k_proxy = len(p_list)
    depths = arange( d_min, d_max+d_by, step=d_by)
    
    solns = zeros((T,K))
    
    for t in range(T):
        solns[t,:] = out[t,0] + cumsum( d_by*out[t,1:(K+1)])
    
    CalA = arange( ceil(solns.min()), floor(solns.max()), step=y_by)
    solns_p = zeros((T,CalA.size))
    
    fig, ax = subplots( nrows=k_proxy+1, ncols=1, sharex=True) #, sharey=True)
    PlotSolnsMC( depths, solns, ax=ax[0])
    ax[-1].set_xlabel("Depth (cm)")
    ax[0].set_ylabel("Age (y BP)")
    
    fig_p, ax_p = subplots( nrows=k_proxy, ncols=1, sharex=True)
    ax_p[-1].set_xlabel("Age (y BP)")
    
    rt = [CalA]
    for ax_n,p in enumerate(p_list):
        inter_p = interp1d( proxy[:,0], proxy[:,p]) 
        print("\nProxy (%s)" % (proxy_names[p]))
        for i,g in enumerate(CalA):
            print(g, end=' ')
            for t in range(T):
                tmp = where(solns[t,:] < g)[0]
                if len(tmp) == 0:
                    indx = K-2
                else:
                    indx = tmp[-1]
                d = depths[indx] + (g-solns[t,indx])/out[t,indx+2]
                if (proxy[0,0] <= d) and (d <= proxy[-1,0] ):
                    solns_p[t,i] = inter_p(d)
                else:
                    solns_p[t,i] = None
    
        ax[ax_n+1].plot( proxy[:,0], proxy[:,p], '-')
        ax[ax_n+1].set_ylabel("Proxy (%s)" % (proxy_names[p]))
        
        PlotSolnsMC( CalA, solns_p, ax=ax_p[ax_n])
        ax_p[ax_n].set_ylabel("Proxy (%s)" % (proxy_names[p]))
        rt += [solns_p]
    
    fig.tight_layout()
    fig_p.tight_layout()
    return rt


if __name__ == "__main__":
    ProxyUQ( "HP1C/HP1C", p_list=[])
    #ProxyUQ( core, p_list=[2,3], y_by=1)
    ProxyUQ( "Auassat/Auassat", p_list=[])
    CalA, solns_AC, solns_CAR = ProxyUQ( "Auassat/Auassat", p_list=[2,3], y_by=10)
    #PlotSolnsMC( CalA, solns_AC)
    fig, ax = subplots()
    ax.hist(solns_AC[25,:]) #Proxy at age CalA[25]
    ax.set_xlabel("Acc. Carbon")
    ax.set_title("%.0f, y BP" % (CalA[25]))

