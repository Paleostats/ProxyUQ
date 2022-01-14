#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 09:41:10 2022

@author: jac
"""

from numpy import loadtxt, arange, zeros, cumsum, where, ceil, floor, array, ndarray
from scipy.interpolate import interp1d

from matplotlib.pylab import plot, subplots, hist
from pandas import read_csv
import json

from pyplotfuncmc import PlotSolnsMC

class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, ndarray):
                return obj.tolist()
            return json.JSONEncoder.default(self, obj)

def ProxyUQ(core, folder="Cores", p_list=None, y_by=10, saveFile = None, plotProxy = None):
    """Produce a posterior MC sample of a proxy, from the MC output of either Plum or Bacon
       where the output is a properly formatted `.out` file. 
       core -- Core handle, e.g. "LL14", associated with the set of `_settings.txt`,
       folder -- Folder for all core output files.  Defaults to "Cores".
       `_proxy.txt`, and `.out` files of interest.  These files should be located within
       the folder identified using the parameter `folder`.
       p_list -- An array of integers identifying the columns within the `core_proxy.txt` file
       to be used to produce the posterior estimates.  If `p_list is None` (default) then 
       the function return the column headings in the `_proxy.txt` file.
       y_by -- The date interval along which estimates will be interpolated, in years.
       returns a `dict` object with three elements: `grid`, the age grid,
                 and an MC sample for each proxy analyzed in `proxy`.
       saveFile -- string of the filename for dict object.
       plotProxy -- string of filename in which plots will be stored. 
    """
    
    ### Read the first three lines of settings
    # Check that the file exists:
    try:
        d_min, d_max, d_by = loadtxt("%s/%s_settings.txt" % (folder, core), max_rows=3)
    except OSError as ex:
        print("The file %s/%s_settings.txt does not exist in this location." % (folder, core))
        return None
    
    K1 = int(floor((d_max-d_min)/d_by))+1
    
    try:
        out = loadtxt( "%s/%s_%d.out" % (folder, core, K1), delimiter=' ')
    except OSError as ex:
        print("The file %s/%s_%d.out does not exist in this location." % (folder, core, K1))
        return None
    
    T, tmp = out.shape
    K = tmp - 3
    
    ### Read the proxie file and, for the moment, change missing values to 0
    try:
        proxy = read_csv("%s/%s_proxies.csv" % (folder, core))
    except OSError as ex:
        print("The file %s/%s_proxies.csv does not exist, do you have a proxy file?" % (folder, core))
        return None
    
    proxy = proxy.fillna(0)
    proxy_names = list(proxy.columns)
    proxy = array(proxy)
    
    if p_list is None:
        print("\n%s_%d.out" % (core,K1))
        print("d_min=%f, d_max=%f, d_by=%f, K=%d, T=%d" % (d_min,d_max,d_by,K,T))
        print("Proxies: ", proxy_names, "\n")
        return None
    
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
    rt = {'dates': CalA}
    proxyResult = []
    
    ### Iterate through the proxy list
    for ax_n,p in enumerate(p_list):
        ### Interpolate the proxy
        inter_p = interp1d( proxy[:,0], proxy[:,p])
        #print("\nProxy (%s)" % (proxy_names[p]))
        ### Itearete through the calendar ages
        for i,g in enumerate(CalA):
            ### For each year, iterate through all the age-depth models
            #print(g, end=' ')
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
        #print('\n')
        ax[ax_n+1].plot( proxy[:,0], proxy[:,p], '-')
        ax[ax_n+1].set_ylabel("Proxy (%s)" % (proxy_names[p]))
        ### Plot the proxy ate cal age with quentiles       
        PlotSolnsMC( CalA, solns_p, ax=ax_p[ax_n,0])
        ax_p[ax_n,0].set_ylabel("Proxy (%s)" % (proxy_names[p]))
        proxyResult += [solns_p.copy()]
    
    rt['proxy'] = {}

    for j in range(len(p_list)):
        rt['proxy'][proxy_names[p_list[j]]] = proxyResult[j]
        
    if saveFile is not None:
        try:
            with open(saveFile, 'w') as outfile:
                json.dump(rt, outfile, cls = NumpyEncoder)
        except Exception as ex:
            print(ex)
            return None

    if plotProxy is not None:
        try:
            fig.tight_layout()
            filename = plotProxy+".png"
            fig.savefig(filename)
            fig_p.tight_layout()
            filename = plotProxy+"_p.png"
            fig_p.savefig(filename)
        except Exception as ex:
            print(ex)
            return None

    return rt

if __name__ == "__main__":
       ### Read the files and provide general info    aa=ProxyUQ( "HP1C/HP1C", p_list=[1], y_by=1)
#    ProxyUQ( "Auassat/Auassat", p_list=[])
    
    aa = ProxyUQ( core = "LL14", folder = "AtmosphericRivers", p_list=[1], saveFile="myjson.json")
    
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
