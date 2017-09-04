''' mod_nc
    Author: Alex Paterson
    
Contains functions for postprocessing of the GFDL output .nc datafiles
'''
from __future__ import division
import numpy as np
from useful import *

import os
import sys


def tropopause_calc(ds_mean):
    '''run the zonal mean tropopause height calculation'''
    i = 0
    h_trop_calc = np.zeros((ds_mean.pentad.size, ds_mean.lat.size))
    if 'pentad' in ds_mean.keys():
        mean_temp = ds_mean.temp.isel(pentad=-1).mean().values
        for k in range(ds_mean.pentad.size):
            for j in range(ds_mean.lat.size):
                h_trop_calc[k,j] = calc_tropopause(ds_mean.temp.isel(pentad=k,lat=j).values, ds_mean.pfull.values, mean_temp)
                i += 1
                percent = i/h_trop_calc.size*100
                sys.stdout.write("\rTropopause: %.2f%%" % percent)
                sys.stdout.flush()            
    else:
        mean_temp = ds_mean.temp.isel(time=-1).mean().values
        for k in range(ds_mean.time.size):
            for j in range(ds_mean.lat.size):
                h_trop_calc[k,j] = calc_tropopause(ds_mean.temp.isel(time=k,lat=j).values, ds_mean.pfull.values, mean_temp)
                i += 1
                percent = i/h_trop_calc.size*100
                sys.stdout.write("\rTropopause: %.2f%%" % percent)
                sys.stdout.flush()
    ds_mean['h_trop_calc'] = (ds_mean.ps.dims, h_trop_calc/1000)
    ds_mean['p_trop_calc'] = (ds_mean.ps.dims, 1000*np.exp(-9.81/(287.04*mean_temp)*h_trop_calc))
    ds_mean.coords['h'] = 287.04*mean_temp/9.81*np.log(1000/ds_mean.pfull)
    return ds_mean


def add_times(ds):
    ds.time.data = ds.time.data/24 - 1
    ds.coords['day'] = ds.time%360
    lag = calc_seasonal_lag(ds.teq.mean('lon'))
    if np.isnan(lag):
        lag = 0   
    ds.coords['day0'] = (ds.day-lag+45)%360
    ds.coords['season'] = np.floor(ds.day0/90)
    ds.coords['pentad'] = np.floor(ds.day0/5)
    return ds

def main_climatology():
    data_dir = '/scratch/ap587/dry_data/'
    exp = 'exp5'
    exps = []
    for name in os.listdir(data_dir):
        if name.split('_')[0] == exp:
            exps += [name]
     
    exps = ['exp9v1_hc10.0-obl0.0', 'exp9v1_hc10.0-obl15.0', 'exp9v1_hc10.0-obl30.0', 'exp9v1_hc10.0-obl45.0', 'exp9v1_hc10.0-obl60.0', 'exp9v1_hc10.0-obl75.0', 'exp9v1_hc10.0-obl90.0']
    overwrite = True
    
    save_dir = '/scratch/ap587/dry_data/climatologies/'
    ncfile='climatology.nc'
    runs = range(0,49)
    for exp in exps:
        if (not overwrite) and (os.path.isfile(save_dir+exp+'/'+ncfile)):
            sys.stdout.write("Skipping: file " + exp + ' exists and overwrite is False\n')            
        else: 
            sys.stdout.write("Processing climatology for experiment " + exp + '\n')
            ds = open_runset(data_dir, exp+'/', runs)
            ds = ds.isel(time=slice(720,None))
            ds = add_times(ds)        
            ds = ds.groupby('pentad').mean('time')
            sys.stdout.write("Dataset loaded, computing EMF\n")
            ds = calc_emf(ds)
            ds_mean = ds.mean('lon')
            ds_mean.load()
            sys.stdout.write("Now computing mass streamfunction\n")
            if ('h_trop' in ds.keys()):
                ds_mean['p_trop'] = 1000*np.exp(-9.81/(287.04*ds_mean.temp.isel(pentad=-1).mean().values)*ds_mean.h_trop*1000) 
                ds_mean = tropopause_calc(ds_mean)
            ds_mean = psi_calc(ds_mean)
            
            save_runset(save_dir, exp, ds_mean, ncfile)
                    


def main():
    data_dir = '/scratch/ap587/dry_data/'
    exp = 'exp9v1'
    exps = []
    for name in os.listdir(data_dir):
        if name.split('_')[0] == exp:
            exps += [name]
    # old code for running specific parts of an experiment
    #vals = [0.1, 0.5, 1, 2, 4, 8]
    #exps = [exp+'%.1f' % val for val in vals]
    #exps = ['exp4_HS_hc2000.0']
    
    save_dir = '/scratch/ap587/dry_data/processed/'
    ncfile='daily_zmean.nc'
    runs = range(0,49)
    for exp in exps:
        sys.stdout.write("Processing experiment " + exp + '\n')
        ds = open_runset(data_dir, exp+'/', runs)
        #ds = emf_calc(ds)
        ds_mean = ds.mean('lon')
        ds_mean.load()
        sys.stdout.write("Dataset loaded, computing\n")
        if ('h_trop' in ds.keys()):
            ds_mean['p_trop'] = 1000*np.exp(-9.81/(287.04*ds_mean.temp.isel(time=-1).mean().values)*ds_mean.h_trop*1000) 
            ds_mean = tropopause_calc(ds_mean)
        ds_mean = psi_calc(ds_mean)        
        ds_mean.time.data = ds_mean.time.data/24 - 1
        ds_mean.coords['year'] = ds_mean.time//360
        ds_mean.coords['day'] = ds_mean.time%360
        lag = calc_seasonal_lag(ds_mean.teq)
        if np.isnan(lag):
            lag = 0
        ds_mean.coords['day_cor'] = (ds_mean.day-lag+45)%360
        ds_mean.coords['season'] = np.floor(ds_mean.day_cor/90)
        
        save_runset(save_dir, exp, ds_mean, ncfile)

        
if __name__ == "__main__":
    main_climatology()