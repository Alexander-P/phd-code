
from __future__ import division
import numpy as np
import scipy.ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import dask
import os.path
matplotlib.rcParams['figure.figsize'] = (12.0, 8.0)
matplotlib.rcParams['contour.negative_linestyle']= 'dashed'
plt.rcParams['image.cmap'] = 'RdBu'

from progress_bar import *
from useful import *

global data_dir
data_dir = '/scratch/ap587/dry_data/'




# --- Attempting to produce plots similar to Ruth's of log(ITCZ lat) v log(500hpa psi)
#     over 5d means for some experiment ---

plot_dir = '/scratch/ap587/plots/ITCZvPSIv2/'
#exp = 'exp5_obl60.0'

exp = 'exp9'
exps = get_exps(data_dir+'processed/',exp)

for exp in exps:
    print('Doing {}'.format(exp))
    ds = open_zmean_runset(data_dir,exp)
    ds = ds.isel(time=slice(720,None))
    
    lag = calc_seasonal_lag(ds.teq)
    if np.isnan(lag):
        lag = 0

    ds.coords['day'] = (ds.day-lag+45)%360
    ds.coords['season'] = np.floor(ds.day/90)
    ds.coords['d5'] = np.floor(ds.day/5)

    ds = ds.groupby('d5').mean('time')

    itcz = np.zeros(len(ds.d5))
    strength = np.zeros(len(ds.d5))

    for i in range(len(ds.d5)):
        data = ds.isel(d5=i)
        #print('5d %d' %i)    
        #edges = get_hc_edges(data)
        #print(edges)
        threshold = 0.15*np.max(np.abs(data.psi.sel(pfull=500,method='nearest').values))
        edges = get_hc_edges(data, method='psi',threshold=threshold)
    
        #print(edges)
        #print('-----')
        itcz[i] = np.abs(edges[1])
        itcz[i] = np.abs(calc_itcz_lat(data))
        strength[i] = np.abs(data.psi).sel(pfull=500,method='nearest').max().values
    
    
    #plot stuff is good
    eq = np.append(np.arange(0,18), np.arange(36,54))
    sw = np.append(np.arange(18,36), np.arange(54,72))
    d = np.arange(0,72)
    cax = np.cos((d-8)*2*np.pi/36)
    a = np.repeat(0, 9)
    a = np.append(a, np.repeat(1, 18))
    a = np.append(a, np.repeat(0, 9))
    marks = np.append(a, a)
    
    a = np.repeat(0, 27)
    a = np.append(a, np.repeat(1, 36))
    #marks = np.append(a, np.repeat(0, 9))
    
    fig, ax = plt.subplots(1,1)
    #ax.plot(itcz[eq], strength[eq], 'bo', label='Equinox')
    #ax.plot(itcz[sw], strength[sw], 'ro', label='Winter')
    f = ax.scatter(itcz[marks==0], strength[marks==0], c=cax[marks==0],
               marker='o', s=40)#, label='Equatorward')
    f2 = ax.scatter(itcz[marks==1], strength[marks==1], c=cax[marks==1],
               marker='^', s=40)#, label='Poleward')
    
    c1 = ax.plot((0,1),(0,0), color='b', marker='s', linestyle='', label='Spring & Autumn')
    c2 = ax.plot((0,1),(0,0), color='r', marker='s', linestyle='', label='Summer & Winter')
    c1 = ax.plot((0,1),(0,0), color='k', marker='o', linestyle='', label='Hotspot equatorward')
    c2 = ax.plot((0,1),(0,0), color='k', marker='^', linestyle='', label='Hotspot poleward')

    #ax.set_xscale('log')
    ax.set_xlim([np.min(itcz)-0.1, np.max(itcz)+0.5*np.max(itcz)])
    ax.set_xticks([0.625, 1.25, 2.5, 5, 10, 20, 40])
    ax.set_xticks(np.arange(0,41,10))
    ax.set_xlim([0,40])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ax.set_yscale('log')
    ax.set_ylim([np.min(strength)-5e9, np.max(strength)+0.3*np.max(strength)])
    ax.set_yticks(np.linspace(np.ceil(ax.get_ylim()[0]/1e10)*1e10, np.floor(ax.get_ylim()[1]/1e10)*1e10, 3))
    #ax.set_yticks([2e10, 4e10, 8e10, 16e10])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ax.set_title(exp[5:], fontsize=14)
    ax.set_xlabel('Ascending branch latitude')
    ax.set_ylabel('PSI_500')
    fig.colorbar(f)
    plt.legend(loc=2, title='', fontsize=8)
    plt.grid()
    
    figname = exp + '_ITCZvPSI'
    plt.savefig(plot_dir + figname + '.pdf', format='pdf')
    plt.savefig(plot_dir + figname + '.png', format='png')
    plt.close()