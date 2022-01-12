#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 11:33:08 2021

@author: jac
"""

from numpy import zeros, quantile, append, array, isnan
#from scipy.interpolate import interp1d
from matplotlib.pylab import plot, hist, subplots


def PlotFuncAtMC( f, x, sample, ax=None, color='blue', histtype='bar'):
    """Plot a histogram of $f_\theta(x)$
        for every $\theta$ in MC sample.
    """
    if ax == None:
        fig, ax = subplots()
    
    T = sample.shape[0] # sample size
    solns = zeros(T)
    for i in range(T):
        solns[i] = f( x, pars=sample[i,:])
    ax.hist( solns, color=color,histtype=histtype )
    return ax, solns


def PlotFuncEvolveMC( f, x, sample, q=[0.10,0.25], ax=None, fill=True, color='blue', med_col='red'):
    """Plot evolution of quantiles q of $f_\theta(x)$
        for every $\theta$ in MC sample.
        q: quantiles below 0.5, these will be repeated with 1-q to make quantiles ranges.
        fill: If True, fill the quantile rages with color, otherwiaw, draw color lines for each quantile.
        med_col: color to plot the median.
    """

    if ax == None:
        fig, ax = subplots()
    
    q = array(q)
    q = append(append( q, [0.5]), 1-q)
    ### First, evaluate the function
    T = sample.shape[0] # sample size
    m = x.size #Number of points in x
    solns = zeros((T,x.size))
    for i in range(T):
        solns[i,:] = f( x, pars=sample[i,:])
    n = q.size
    quan = zeros((n,m))
    ### Calculate the quantiles
    for j in range(m):
        quan[:,j] = quantile( solns[:,j], q)
    ###Plot the quiantiles:
    k = n//2
    if fill:
        for i in range(k):
            ax.fill_between( x, quan[i,:], quan[-1-i,:], color=color, alpha=0.5/k)
    else:
        for i in range(n):
            ax.plot( x, quan[i,:], '--', color=color, linewidth=0.5)
    ###Plot the median
    ax.plot( x, quan[k,:], '-', color=med_col, linewidth=1.5)
    
    return ax, quan


def PlotSolnsMC( x, solns, q=[0.10,0.25], ax=None, fill=True, color='blue', med_col='red'):
    """Plot evolution of quantiles q of $f_\theta(x)$
        for every $\theta$ in MC sample.
        q: quantiles below 0.5, these will be repeated with 1-q to make quantiles ranges.
        fill: If True, fill the quantile rages with color, otherwiaw, draw color lines for each quantile.
        med_col: color to plot the median.
    """

    if ax == None:
        fig, ax = subplots()
    
    q = array(q)
    q = append(append( q, [0.5]), 1-q)
    n = q.size
    
    m = x.size #Number of points in x
    T, m1 = solns.shape #sample size and number of points in x
    if not(m == m1):
        print("Size of x (%d) not equal to the number of columns (%d) in solns" %\
              (m,m1))
        return
    
    quan = zeros((n,m))
    ### Calculate the quantiles
    for j in range(m):
        tmp = solns[:,j]
        tmp = tmp[~isnan(tmp)]
        if tmp.size > 100:
            quan[:,j] = quantile( tmp, q)
        else:
            quan[:,j] = [None]*n
    ###Plot the quiantiles:
    k = n//2
    if fill:
        for i in range(k):
            ax.fill_between( x, quan[i,:], quan[-1-i,:], color=color, alpha=0.5/k)
    else:
        for i in range(n):
            ax.plot( x, quan[i,:], '--', color=color, linewidth=0.5)
    ###Plot the median
    ax.plot( x, quan[k,:], '-', color=med_col, linewidth=1.5)
    
    return ax, quan

