#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 19:57:31 2021

@author: rosswilkinson

SUMMARY: 
    Creates a violin plot of a 1-d array of data.
INPUT:
    meas - 1-d numpy array.
OUTPUT:
    No direct output. Creates new figure or is added to the current axes.
"""
import numpy as np
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt

def violinPlot(meas, marker, ax):
    
    ### Create probability distribution of data
    mu = meas.mean()
    sigma = meas.std()
    x = np.linspace(mu-(sigma*4), mu+(sigma*4),1000)
    
    # Guassian kernel density estimate
    gkde = stats.gaussian_kde(meas)
    # Evaluate at x    
    kdepdf = gkde.evaluate(x)
    
    ### Calculate 95% CI
    n = len(meas)
    df = n - 1
    tstar = stats.t.ppf([0.025,0.975],df)
    se = sigma / np.sqrt(n)
    ci95 = mu + (tstar*se)
    es = mu / sigma
    
    ### Fill area of vertical mirrored distribution
    c = np.array([1,1,1])*0.8
    ax.fill_betweenx(x, kdepdf, -kdepdf, color=c)
    
    for data in meas:
        ax.scatter(0,data,marker = marker)
    
    kLow = np.argmax(x>ci95[0])
    kHigh = np.argmax(x>ci95[1])
    kMu = np.argmax(x>mu)
    
    ax.plot([-kdepdf[kLow], kdepdf[kLow]], [x[kLow], x[kLow]], linestyle = ':', linewidth = '1', color = 'k')
    ax.plot([-kdepdf[kHigh], kdepdf[kHigh]], [x[kHigh], x[kHigh]], linestyle = ':', linewidth = '1', color = 'k')
    ax.plot([-kdepdf[kMu], kdepdf[kMu]], [x[kMu], x[kMu]], linestyle = '-', linewidth = '1', color = 'k')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    return kdepdf, x, ci95, mu, sigma, es