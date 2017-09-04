# useful functions for data analysis
from __future__ import division
import numpy as np
import os.path
import xarray as xr
import matplotlib.pyplot as plt
import gradients as gr


def file_exists(filename):
    import os.path
    return os.path.isfile(filename)

def rfile(data_dir, exp, run, ncfile='daily.nc'):
    """Returns the file for a given ratio and month"""
    import os.path
    return os.path.join(data_dir, exp, 'run%d' % run, ncfile)

def rfile00(data_dir, exp, run, ncfile='daily.nc'):
    """Returns the file for a given ratio and month"""
    import os.path
    return os.path.join(data_dir, exp, 'run%.3d' % run, ncfile)

def open_runset(data_dir, exp, runs=range(100)):
    fnames = [rfile(data_dir, exp, d) for d in runs] + [rfile00(data_dir, exp, d) for d in runs]
    fnames = [f for f in fnames if file_exists(f)]
    if fnames:
        import xarray as xr
        return xr.open_mfdataset(fnames, chunks={'time': 90}, decode_times=False)
    else:
        print('WARNING: ' + data_dir + exp + ' not found')
        print([rfile00(data_dir, exp, d) for d in runs])
        return None
    
def open_zmean_runset(data_dir, exp):
    import os.path
    fnames = os.path.join(data_dir, 'processed', exp, 'daily_zmean.nc')
    if fnames:
        import xarray as xr
        return xr.open_mfdataset(fnames, chunks={'time': 90}, decode_times=False)
    else:
        return None
    
def open_climatology(data_dir, exp):
    import os.path
    fnames = os.path.join(data_dir, 'climatologies', exp, 'climatology.nc')
    if fnames:
        import xarray as xr
        return xr.open_mfdataset(fnames, decode_times=False)
    else:
        return None
    
def save_runset(save_dir, exp, ds, ncfile='daily.nc'):
    import os.path
    dest = os.path.join(save_dir, exp + '/')
    if not os.path.exists(dest):
        os.makedirs(dest)
    if dest:
        ds.to_netcdf(dest + ncfile)
        
def get_exps(data_dir, exp):
    #return a list of all exps in the data_dir folder with the argument as the first section
    exps = []
    vals = []
    for name in os.listdir(data_dir):
        if name.split('_')[0] == exp:        
            exps += [name]
            vals += [name.split('_')[1]]
    return exps

def get_vals(data_dir, exp):
    #return a list of all experimental values in the data_dir folder with the argument as the first section
    exps = []
    vals = []
    for name in os.listdir(data_dir):
        if name.split('_')[0] == exp:        
            exps += [name]
            vals += [name.split('_')[1]]
    return vals
    
    
def plot_grid(data, title_text='', nlev=15, vmint=120, vmaxt=300, vmaxu=50, vmaxv=1.2, vmaxsf=0.5e11, do_log=0, do_tp=1):
    '''Inputs: data -- 3d grid with axes lon, lat, pfull
    Plot produced: subplots of 4 quantities'''
    
    h_trop_exists = 'p_trop' in data.keys()
    tp_calc_exists = 'p_trop_calc' in data.keys()
    psi_exists = 'psi' in data.keys()
    
    if do_log:
        data['pfull'] = np.log(data.pfull)
        if h_trop_exists:
            data['p_trop'] = np.log(data.p_trop)
        if tp_calc_exists:
            data['p_trop_calc'] = np.log(data.p_trop_calc)
    
    fig, axes = plt.subplots(ncols=2,nrows=2,figsize=(21,12),sharey=True)
    if psi_exists:
        data.teq.plot.contourf(ax=axes[0,0], levels=nlev,vmax=vmaxt,vmin=vmint, extend='both')
        data.temp.plot.contourf(ax=axes[0,1],levels=nlev,vmax=vmaxt,vmin=vmint, extend='both')
        data.ucomp.plot.contourf(ax=axes[1,0],levels=nlev,vmax=vmaxu, extend='both')
        data.psi.plot.contourf(ax=axes[1,1],levels=nlev,vmax=vmaxsf, extend='both')
    else:
        data.teq.plot.contourf(ax=axes[0,0], levels=nlev,vmax=vmaxt,vmin=vmint, extend='both')
        data.temp.plot.contourf(ax=axes[0,1],levels=nlev,vmax=vmaxt,vmin=vmint, extend='both')
        data.ucomp.plot.contourf(ax=axes[1,0],levels=nlev,vmax=vmaxu, extend='both')
        data.vcomp.plot.contourf(ax=axes[1,1],levels=nlev,vmax=vmaxv, extend='both')    
    
    if (h_trop_exists):
        data.p_trop.plot.line('k',ax=axes[0,0])
        data.p_trop.plot.line('k',ax=axes[0,1])
        data.p_trop.plot.line('k',ax=axes[1,0])
        data.p_trop.plot.line('k',ax=axes[1,1])   
        
    if (tp_calc_exists and do_tp):
        data.p_trop_calc.plot.line('k--',ax=axes[0,0])
        data.p_trop_calc.plot.line('k--',ax=axes[0,1])
        data.p_trop_calc.plot.line('k--',ax=axes[1,0])
        data.p_trop_calc.plot.line('k--',ax=axes[1,1])  
  
    axes[0,0].set_title('Equilibrium temperature')
    axes[0,1].set_title('Model temperature')
    axes[1,0].set_title('Zonal wind')
    if psi_exists:
        axes[1,1].set_title('Mass streamfunction')
    else:
        axes[1,1].set_title('Meridional Wind')
    fig.suptitle(title_text, fontsize=20)    
    plt.gca().invert_yaxis()

def plot_timeseries(ds, title_text='', nlev=15, vmin=150, vmax=300):
    fig, axes = plt.subplots()
    ds.isel(pfull=-1).plot.contourf(x='day', y='lat', ax=axes, levels=nlev, vmin=vmin, vmax=vmax, extend='both')
    plt.title(title_text)
    
    
def plt_vert_structure(runset):
    fig, axes = plt.subplots(nrows=2,ncols=2, sharey=True, sharex=True)    
   
    sel_data = runset.mean('lon').mean('time').isel(lat=32)
    #vert_coord = sel_data.pfull.data
    vert_coord = 287.04*255/9.81*np.log(1000/sel_data.pfull.data)
    axes[0,0].plot(sel_data.teq.data, vert_coord, 'r--', label = 'Equatorial teq')
    axes[0,0].plot(sel_data.temp.data, vert_coord, 'r-', label = 'Equatorial temp')
    
    sel_data = runset.mean('lon').mean('time').isel(lat=21)
    axes[0,1].plot(sel_data.teq.data, vert_coord, 'b--', label = 'Subtropical teq')
    axes[0,1].plot(sel_data.temp.data, vert_coord, 'b-', label = 'Subtropical temp')

    sel_data = runset.mean('lon').mean('time').isel(lat=10)
    axes[1,0].plot(sel_data.teq.data, vert_coord, 'k--', label = 'Midlatitudes (60N) teq')
    axes[1,0].plot(sel_data.temp.data, vert_coord, 'k-', label = 'Midlatitudes (60N) temp')
    
    sel_data = runset.mean('lon').mean('time').isel(lat=1)
    axes[1,1].plot(sel_data.teq.data, vert_coord, 'g--', label = 'Polar (85N) teq')
    axes[1,1].plot(sel_data.temp.data, vert_coord, 'g-', label = 'Polar (85N) temp')
    
    axes[1,1].set_xlabel('Temperature (K)')
    axes[1,0].set_xlabel('Temperature (K)')
    axes[0,0].set_ylabel('Pressure (hPa)')
    axes[1,0].set_ylabel('Pressure (hPa)')
    
    axes[0,0].set_title('Tropics')
    axes[0,1].set_title('Subtropics')
    axes[1,0].set_title('Midlatitudes')
    axes[1,1].set_title('Polar')
    
    #plt.gca().invert_yaxis()
    plt.tight_layout()
    
    
def calc_seasonal_lag(data, method='cross_eq'):
    """Return the index of the first equinox - where temperature max is at the equator - over the input time series.
    or where the temperature max crosses the equator"""
    surf_t = data.isel(pfull=-1)
    
    if method == 'at_eq':
        max_lati = np.argmax(surf_t.isel(time=0).data)
        i = 0
        while ((max_lati != (len(surf_t.lat)//2)) and (i < 360)):
            i += 1
            max_lati = np.argmax(surf_t.isel(time=i).data)
        if i==360:
            print('WARNING: No seasonal lag found, returning nan')
            lag = np.nan
        else:
            lag = i
        
    elif method == 'cross_eq':
        lati = np.argmax(surf_t.isel(time=0).data)
        latm1 = data.lat.isel(lat=lati)
        exit = False
        i = 0
        while (not exit) and (i<360):
            i += 1
            lati = np.argmax(surf_t.isel(time=i).data)
            lat = data.lat.isel(lat=lati)
            if np.sign(lat.data) != np.sign(latm1.data):
                exit = True
                lag = i
            elif i==360:
                print('WARNING: No seasonal lag found, returning nan')
                lag = np.nan
            latm1 = lat
    
    return lag

def calc_high_summer(data):
    """Return the index of the data field where the maximum data value is furthest North
    Inputs: data - xarray of global zonal mean temperature data"""
    llm = data.isel(pfull=slice(-6,None)).mean('pfull').groupby('day').mean('time')
    peak_lats = np.zeros(llm.day.size)
    for i in range(llm.day.size):
        peak_lats[i] = np.argmax(llm.isel(day=i))
    a, count_no = np.unique(peak_lats, return_counts=True)
    return int(np.argmax(peak_lats) + round(count_no[-1]/2))

def calc_tropopause(temp, p, gm_t=270):
    """Calculate an analytical tropopause height from 1d temperature and pressure
    gm_t is a global mean temperature used for the height-pressure conversion
    Use method descibed in Reichler et al. (2003)"""
    h_trop = np.nan
    crit_lapse = 2/1000     # critical lapse rate (K/m)
    r_depth = 2000          # required depth (m) if lapse rate is below crit_lapse for this depth, tropopause is defined
    k = 2/7
    g = 9.81
    rdgas = 287.05        
    phalfk = np.zeros(len(p)-1)
    heights = np.zeros(len(p)-1)
    lhalf = np.zeros(len(p)-1)
    for i in range(len(p)-1): 
        phalfk[i] = (p[i]**k + p[i+1]**k)/2
        heights[i] = rdgas*gm_t/g*np.log(1000/phalfk[i]**(1/k))
        lhalf[i] = ((temp[i+1]-temp[i])*(p[i]**k+p[i+1]**k)) / ((p[i+1]**k-p[i]**k)*(temp[i]+temp[i+1])) * k*g/rdgas    
    j = 0
    no_trop = 1
    #reverse array dir to iterate from bottom to top
    phalfk = phalfk[::-1]
    lhalf = lhalf[::-1]
    heights = heights[::-1]
    has_been_pos = 0    
    while (no_trop) and (j<len(phalfk)):
        #print('searching for tropopause: %d' % j) 
        #print('lapse rate is: %.4f' % lhalf[j])
        if (lhalf[j] <= crit_lapse) & has_been_pos:        
            is_within_2k = (heights-heights[j]) < r_depth
            is_above =  heights-heights[j] > 0
            ind = is_within_2k & is_above
            is_below_cl = lhalf[ind] < crit_lapse
            if is_below_cl.all():            
                no_trop = 0
                p_tropk = phalfk[j-1] + (phalfk[j] - phalfk[j-1])/(lhalf[j] - lhalf[j-1])*(crit_lapse-lhalf[j-1])
                h_trop = rdgas*gm_t/g*np.log(1000/p_tropk**(1/k))
                #print('Found tropopause, height: %.2f' % h_trop)
        if lhalf[j] > 0:
            has_been_pos = 1
        j += 1
    if h_trop <= 0:
        h_trop = np.nan
    if h_trop == np.nan:
        print('Warning: no tropopause found, returning NaN')
    return h_trop

def check_spinup(data):
    '''look at the temperature yearly means to check for equilibrium'''
    return

def calc_pot_temp(ds):
    p0 = 1000
    R_cp = 0.286    # for dry air
    ds['pot_temp'] = (ds.temp.dims, ds.temp*(p0/ds.pfull)**R_cp)
    return ds

def calc_bv_freq(ds):
    #brunt-vaisala frequency
    if not('pot_temp' in ds.keys()):
        ds = calc_pot_temp(ds)
    if not('height' in ds.keys()):
        ds['height'] = 287.04*ds.temp.mean()/9.81*np.log(1000/ds.pfull)
    dth = ds.pot_temp.diff('pfull')
    dz = ds.height.diff('pfull')
    bv = np.sqrt(9.81/ds.pot_temp*dth/dz)
    ds['bv'] = bv
    return ds

def calc_eady_growth(ds):
    if not('bv' in ds.keys()):
        ds = calc_bv_freq(ds)
    f = 2*7.2921e-5*np.sin(np.deg2rad(ds.lat))
    du = ds.ucomp.diff('pfull')
    dz = ds.height.diff('pfull')
    ds['eady'] = 0.3098*9.81*abs(f)*abs(du/dz)/ds.bv
    return ds

def running_mean(y, dim=0, window=1):
    '''Calculate a mean along dimension dim of y with a window of specified size -- required odd'''
    if window%2 == 0:
        print('window needs to be odd')
    cumsum = np.cumsum(np.insert(y, 0, 0))
    return (cumsum[window:] - cumsum[:-window]) / window

def nan_running_mean(y, dim=0, window=1):
    '''running mean but padded to original shape with np.nan'''

    nan_rm = np.full(y.shape, np.nan)
    rm = running_mean(y, dim, window)    
    nan_rm[window//2:-window//2+1] = rm    
    return nan_rm


def calc_itcz_lat(ds, do_plot=False):
    '''Locate the latitude of the ITCZ using ds.omega, the vertical velocity
        ds: dataarray with dimensions (pfull, lat)'''
    x = ds.lat.values
    y = (-ds.omega).isel(pfull=slice(-16,-1)).sum('pfull').values
   
    run_y = nan_running_mean(y, window=9)
    drun_y = np.gradient(run_y)
    ddrun_y = np.gradient(drun_y)
    roots = calc_roots(x, drun_y)
    
    vals = []
    for root in roots:
        i = (np.abs(x-root)).argmin()
        vals.append((root, i, float(run_y[i]), float(ddrun_y[i])))

    # vals is a list of [root, root_index, omega, dd_omega] lists
    vals.sort(key=lambda tup: tup[2], reverse=True)
    #print vals

    itcz = 90
    for val in vals:
        if (val[3] < 0) and (np.abs(val[0]) < np.abs(itcz)) and (val[2] > 0.5*vals[0][2]):
            itcz = val[0]
    if itcz == 90:
        itcz = np.nan
        
    if (do_plot == True):
        fig,ax = plt.subplots()
        plt.plot(x, y, label='omega')
        plt.plot(x, run_y, label='running mean')
        plt.plot(x, drun_y, label='d(runing)')
        plt.plot(x, ddrun_y, label='dd(runing)')
        [ax.axvline(extreme, color='k', linewidth=0.5) for extreme in roots]
        ax.axvline(itcz, color='r', linewidth=1)
        #ax.set_xlim([-30,30])
        #ax.set_ylim([-0.5,0.5])
        plt.legend(loc=1)
        plt.grid()
        
    if (np.isnan(itcz)):
        print('No ITCZ found, returning nan')
    return itcz


def calc_descending_branch(ds, do_plot=False):
    '''Locate the latitude of the descending branch of the Hadley cell using 
    ds.omega, the vertical velocity.
    Can def merge this with itcz latitude
        ds: dataarray with dimensions (pfull, lat)'''
    x = ds.lat.values
    y = (-ds.omega).isel(pfull=slice(-16,-1)).sum('pfull').values
   
    run_y = nan_running_mean(y, window=9)
    drun_y = np.gradient(run_y)
    ddrun_y = np.gradient(drun_y)
    roots = calc_roots(x, drun_y)
    
    vals = []
    for root in roots:
        i = (np.abs(x-root)).argmin()
        vals.append((root, i, float(run_y[i]), float(ddrun_y[i])))

    # vals is a list of [root, root_index, omega, dd_omega] lists
    vals.sort(key=lambda tup: tup[2], reverse=False)
    #print vals

    itcz = 90
    for val in vals:
        if (val[3] > 0) and (np.abs(val[0]) < np.abs(itcz)) and (val[2] < 0.5*vals[0][2]):
            itcz = val[0]
    if itcz == 90:
        itcz = np.nan
        
    if (do_plot == True):
        fig,ax = plt.subplots()
        plt.plot(x, y, label='omega')
        plt.plot(x, run_y, label='running mean')
        plt.plot(x, drun_y, label='d(runing)')
        plt.plot(x, ddrun_y, label='dd(runing)')
        [ax.axvline(extreme, color='k', linewidth=0.5) for extreme in roots]
        ax.axvline(itcz, color='r', linewidth=1)
        #ax.set_xlim([-30,30])
        #ax.set_ylim([-0.5,0.5])
        plt.legend(loc=1)
        plt.grid()
        
    if (np.isnan(itcz)):
        print('No descending point found, returning nan')
    return itcz
    

######################
'''Hadley cell width calculations'''
######################

def calc_hc_edges(ds, method='psi', threshold=0, do_plot=0):
    '''
    locate the two Hadley cell and return its width
    arguments:
        ds: dataarray with elements 'lat' and another with the same label as method and dimensions (pfull, lat)
        method: string defining which field to use for the HC calculation.
                - 'psi' uses 500hPa meridional mass streamfunction
                - 'ucomp' uses surface zonal winds
        threshold: if set non-zero instead look for roots of abs(y)-threshold
        do_plot: optional, display some diagnostic plots
        
    outputs:
        edges: list of 3 edge latitudes
    '''
    xx = ds.lat.values
    if method in ['psi', 'psi2']:
        if not ('psi' in ds.keys()):
            print('Calculating streamfunction')
            ds = psi_calc(ds.mean('lon'))        
        yy = ds.psi.sel(pfull=500,method='nearest')
    elif method == 'ucomp':
        yy = ds.ucomp.isel(pfull=-1)
    x = np.linspace(min(xx), max(xx), 500)
    y = np.interp(x, xx, yy)
    t = ds.temp.isel(pfull=-1).values
    t_max_lat = ds.lat.isel(lat=np.argmax(t))
    omega = ds.omega.sel(pfull=500,method='nearest').values
    w_max_lat = ds.lat.isel(lat=np.argmin(omega))
    
    if threshold is not 0:
        y = np.abs(y) - threshold
        
    edges = [np.nan, np.nan, np.nan]
           
    if 'psi' in method:
        roots = calc_roots(x,y)
            
        # only include roots with a significantly large peak between them
        r = []
        for i in range(1,len(roots)):
            max_between_roots = np.nanmax(np.abs((yy.where((yy.lat >= roots[i-1]) & (yy.lat <= roots[i])))))
            if max_between_roots > 2*threshold:
                r.append(roots[i-1])
                r.append(roots[i])
        roots = np.unique(r)
                
            
        r = []
        i=0
        while i < len(roots):
            if i == len(roots)-1:
                r += [roots[i]]
            elif np.abs(roots[i] - roots[i+1]) < 10:    # 10 degree minimum cell size
                # require a minima in the middle to combine
                idx = (np.abs(x-roots[i])).argmin()
                if (y[idx] - y[idx+1]) > 0:
                    r += [(roots[i] + roots[i+1])/2]
                    i += 1
                else:
                    r += [roots[i]]
            else:
                r += [roots[i]]
            i += 1
                
        roots = np.array(r)
            
        if method == 'psi':
            # ieq is the index of the root closest to the warmest point of the dataset
            itcz_lat, itcz_i = min((val, ieq) for (ieq, val) in enumerate([abs(root) for root in roots-t_max_lat.values]))
            dpsi = []
            for i in range(len(roots)):
                idx = (np.abs(xx-roots[i])).argmin()
                dpsi.append(yy[idx+1] - yy[idx])
                    
            try:
                if (itcz_i == 0):
                    if (dpsi[0] < 0):
                        edges = [np.nan, roots[0], roots[1]]
                    else:
                        edges = [np.nan, roots[1], roots[0]]
                elif (itcz_i == len(roots)-1):
                    if (dpsi[len(roots)-1] < 0):
                        edges = [roots[-2], roots[-1], np.nan]
                    else:
                        edges = [roots[-1], roots[-2], np.nan]
                else:
                    if (dpsi[itcz_i] < 0):
                        edges = [roots[itcz_i-1], roots[itcz_i], roots[itcz_i+1]]
                    else:
                        print('Hotspot closest root', itcz_i)
                        print('Roots')
                        print(roots)
                        print('PSI gradient')
                        print(dpsi)
            except:
                print('ITCZ index: %d' % itcz_i)
                print(roots)
                raise Exception('Error  finding roots for the Hadley cell')
            
        elif method == 'psi2':
            # use the index of the root closest to the maximum of vertical velocity. Then move towards the hotspot to find ITCZ latitude.
            wpeak_lat, upeak_i = min((val, ieq) for (ieq, val) in enumerate([abs(root) for root in roots-w_max_lat.values]))

            
        
    elif method == 'ucomp':
    
        roots = calc_roots(x, y)
        root, ieq = min((val, ieq) for (ieq, val) in enumerate([abs(root) for root in roots-t_max_lat.values]))
            
        if len(roots) >= 4:
            # we likely have two extra roots within the equatorial zone, need to distinguish them
            signs = calc_roots_sign(x,y)
            im1 = np.nonzero((signs == -1) & (roots < roots[ieq]))
            ip1 = np.nonzero((signs == 1) & (roots > roots[ieq]))
            edges = [roots[im1] , roots[ieq], roots[ip1]]
        elif len(roots) == 3:
            # if the central root is within (obliquity) of the equator
            if abs(roots[ieq]) < 23:
                edges = roots
            else:
                middle = np.ones(len(x)) * max(abs(y))
                middle[int(len(y)/2) - 6: int(len(y)/2.) + 6] = y[int(len(y)/2.) - 6: int(len(y)/2) + 6]
                cen_max = x[np.argmin(abs(middle))]
                other_end = abs(roots-t_max_lat.values).argsort()[1]
                edges = [roots[ieq], roots[other_end], cen_max]
                edges.sort()
        elif len(roots) == 2:
            middle = np.ones(len(x)) * max(abs(y))
            middle[len(y)/2 - 6: len(y)/2 + 6] = y[len(y)/2 - 6: len(y)/2 + 6]
            cen_max = x[np.argmin(abs(middle))]
            other_end = abs(roots-t_max_lat.values).argsort()[1]
            edges = [roots[ieq], roots[other_end], cen_max].sort()
            edges = [roots[ieq-1],roots[ieq],roots[ieq+1]]
        elif len(roots) == 1:
            edges = [np.nan,np.nan,np.nan]
            print("Didn't find enough roots to define the Hadley cell")
        else:
            print i, roots
            raise Exception("Error finding roots for the Hadley cell")
            
    else:
        raise Exception("Error: Invalid Hadley cell method %s" % method)

    return edges
    
    #roots = calc_roots(x,y,do_plot)
    # locate the index of the root closest to the equator
    #root, ieq = min((val, ieq) for (ieq, val) in enumerate([abs(root) for root in roots]))
    #if len(roots) == 2:
    #    width = (abs(roots[1]-roots[0]),np.nan)
    #elif len(roots) > 2:
    #    width = (roots[ieq]-roots[ieq-1], roots[ieq+1]-roots[ieq])
    #else:
    #    width = (np.nan,np.nan)
    #return width

def calc_roots(x,y,do_plot=0):
    '''
    Calculate roots in a data set, avoiding the edges of the domain
    Arguments:
        x: array-like
        y: array-like
    Outputs:
        roots: list of roots
    '''
    
    if len(x) != len(y):
        raise ValueError('x and y must have the same dimensions, currently {} and {}'.format(len(x), len(y)))
    roots = []
    for i, val in enumerate(y[5:-5],start=5):
        if not (np.sign(y[i]) == np.sign(y[i+1])):
            roots += [x[i] + (x[i+1]-x[i])*(y[i]/(y[i]-y[i+1]))]
    if do_plot:
        fig, ax = plt.subplots()
        ax.plot(x,y)
        [ax.axvline(root, color='r') for root in roots]
    return roots
                        
def calc_roots_sign(x,y,do_plot=0):
    '''
    Calculate the sign of roots in a data set, avoiding the edges of the domain
    Arguments:
        x: array-like
        y: array-like
    Outputs:
        signs: list of roots signs
    '''
    roots = []
    signs = []
    for i, val in enumerate(y[5:-5],start=5):
        if not (np.sign(y[i]) == np.sign(y[i+1])):
            signs += [np.sign(y[i+1])]
    return signs


def calc_hc_width(ds, rot=7.3e-5):
    '''
    Held & Hou 1980 scaling for hadley cell width
    Method: take slice over the equatorial region
        static stability is the difference between the ground level and tropopause level potential temp
        del_h is the fractional eq-pole pot temp difference
    '''
    tp_height = ds.h_trop.isel(lat=slice(20,44)).mean('lat')
    tp_pres = ds.p_trop.isel(lat=slice(20,44)).mean('lat')
    ds = calc_pot_temp(ds)
    #static_stab = (ds.pot_temp.isel(lat=slice(20,44)).sel(pfull=tp_pres, method='nearest').mean('lat').mean('time') -                    ds.pot_temp.isel(lat=slice(20,44),pfull=-1).mean('lat').mean('time'))
    p0 = 1000
    R_cp = 0.286    # for dry air
    ds['rbal_pot_temp'] = (ds.teq.dims, ds.teq*(p0/ds.pfull)**R_cp)
    del_h_S = ((ds.pot_temp.isel(lat=slice(30,34)).isel(pfull=-1).mean('lat')-
                   ds.pot_temp.isel(lat=slice(0,4)).isel(pfull=-1).mean('lat'))/
             ds.pot_temp.isel(lat=slice(30,34)).isel(pfull=-1).mean('lat'))
    del_h_N = ((ds.pot_temp.isel(lat=slice(30,34)).isel(pfull=-1).mean('lat')-
                   ds.pot_temp.isel(lat=slice(60,None)).isel(pfull=-1).mean('lat'))/
             ds.pot_temp.isel(lat=slice(30,34)).isel(pfull=-1).mean('lat'))
    del_h = np.max(del_h_S, del_h_N)
    R = (9.81*tp_height*1e3*del_h)/((rot)**2 * (6371e3)**2)
    ccw = (5/3 * R)**0.5
    return np.rad2deg(ccw.values)

def calc_heldhou(ds, rot=7.3e-5):
    
    r_cp = 0.286
    itcz = calc_itcz_lat(ds)
    if not np.isnan(itcz):
        trop_h = ds.h_trop_calc.sel(lat=np.arange(itcz-10,itcz+10,ds.lat[2]-ds.lat[1]), method='nearest').mean('lat').values
        trop_p = ds.p_trop_calc.sel(lat=np.arange(itcz-10,itcz+10,ds.lat[2]-ds.lat[1]), method='nearest').mean('lat').values
        ds = calc_pot_temp(ds)
        ds = ds.isel(pfull=-1)
        del_h_S = ((ds.pot_temp.sel(lat=np.arange(itcz-10,itcz+10,ds.lat[2]-ds.lat[1]), method='nearest').mean('lat')-
                   ds.pot_temp.sel(lat=np.arange(-88,-68,ds.lat[2]-ds.lat[1]), method='nearest').mean('lat'))/
                   ds.pot_temp.sel(lat=np.arange(itcz-10,itcz+10,ds.lat[2]-ds.lat[1]), method='nearest').mean('lat')).values
        del_h_N = ((ds.pot_temp.sel(lat=np.arange(itcz-10,itcz+10,ds.lat[2]-ds.lat[1]), method='nearest').mean('lat')-
                   ds.pot_temp.sel(lat=np.arange(68,88,ds.lat[2]-ds.lat[1]), method='nearest').mean('lat'))/
             ds.pot_temp.sel(lat=np.arange(itcz-10,itcz+10,ds.lat[2]-ds.lat[1]), method='nearest').mean('lat')).values
    
        del_h = np.array((del_h_S, del_h_N))
        R = (9.81*trop_h*1e3*del_h)/((rot)**2 * (6371e3)**2)
        return np.rad2deg((5/3 * R)**0.5)
    else:
        return (np.nan, np.nan)
    
    
    
def peak_emf(ds):
    '''
    Return the latitudes of peak eddy momentum flux, a proxy for Hadley cell width according to some theories
    '''
    emf = ds.emf.sel(pfull=500,method='nearest').mean('time').values
    maxi = np.argmax(emf)
    mini = np.argmin(emf)
    return (ds.lat.isel(lat=mini).values, ds.lat.isel(lat=maxi).values)
    
def psi_calc(ds_mean):
    if 'time' in ds_mean.dims:
        data = {'omega':ds_mean.omega.data, 'vcomp':ds_mean.vcomp.data, 'time':ds_mean.time.data, 'pfull':ds_mean.pfull.data, 'lat':ds_mean.lat.data}
    elif 'pentad' in ds_mean.dims:
        data = data = {'omega':ds_mean.omega.data, 'vcomp':ds_mean.vcomp.data, 'time':ds_mean.pentad.data, 'pfull':ds_mean.pfull.data, 'lat':ds_mean.lat.data}
    else:
        print('No time dimension found in psi_calc. About to break')
        raise Exception('Whoops')
    psi = cal_stream_fn(data)
    ds_mean['psi'] = (ds_mean.omega.dims, psi)
    return ds_mean


def calc_emf(ds):
    '''
    Eddy momentum flux calculation
    Inputs:
        ds: xarray dataset containing at minimum fields (ucomp, vcomp, ucomp_vcomp)
    Outputs:
        ds: modified dataset with emf fields added
    '''
    
    uv = -86400. * (ds.ucomp_vcomp - ds.ucomp*ds.vcomp)
    uv_dy = gr.ddy(uv, uv=True)
    
    uw = -86400. * (ds.ucomp_omega - ds.ucomp*ds.omega)
    uw_dp = gr.ddp(uw)

    
    ds['uv'] = (ds.ucomp.dims, uv)
    ds['uv_dy'] = (ds.ucomp.dims, uv_dy)
    ds['uw'] = (ds.ucomp.dims, uw)
    ds['uw_dp'] = (ds.ucomp.dims, uw_dp)   
    return ds

    
####### functions from penny's github code #########
'''
__author__ = 'Penelope Maher'
'''

def MeanOverDim(data,dim):
    'Assumed data is masked with nans for missing values'
    return np.nansum(data,axis=dim)/np.sum(np.isfinite(data),axis=dim)

#Code Description: 
#   eg see text book: P141-142 of global physical climatology by Hartmann 

#   Psi = ( (2 pi a cos(lat)) /g) integrate ^P _0 [v] dp 
#   [v] = (g/(2 pi a cos (lat))) (partial Phi/ partial p) 
#   [omega] = (-g/(2 pi a^2 cos (lat))) (partial Phi/ partial lat) 

def cal_stream_fn(data):
    'Calculates the mass stream function. \
  Data is a dictionary containing masked arrays of omega, v wind, lat, time and pressure.'

    grav=9.80665      #m/s
    rad_earth=6378100.0     #m
    lat_rad=np.radians(data['lat'])	
    fn_const_dlat = (-2*(np.pi)*(rad_earth*rad_earth))/grav  

    #integrate right to left [omega] and then integrate left to right [omega]
  
    psi_p = np.zeros((len(data['time']),len(data['pfull']),len(data['lat'])))
    psi_lat = np.zeros((len(data['time']),len(data['pfull']),len(data['lat'])))
    psi_p_rev = np.zeros((len(data['time']),len(data['pfull']),len(data['lat'])))
    psi_lat_rev = np.zeros((len(data['time']),len(data['pfull']),len(data['lat'])))
    omega = data['omega']
    vcomp = data['vcomp']


    #assert isinstance(omega,np.ma.core.MaskedArray)  == True , 'Data is not masked' 
  
    print '... Integrate first direction'
  
    for t in xrange(len(data['time'])):
        for p in xrange(len(data['pfull']) -1):  #0 to 34 inclusive
            #print('pressure %d / %d' % (p, len(data['pfull'])))
            psi_lat[t,p,0]=0
            psi_lat_rev[t,p,0]=0
            
            for l in xrange(len( data['lat']) -1):
                w_val = np.array([omega[t,p,l],omega[t,p+1,l],omega[t,p,l+1],omega[t,p+1,l+1]])

                dlat=lat_rad[l] - lat_rad[l+1]   #negative (era in IDL will be opposite order)
                psi_lat[t,p,l+1] = psi_lat[t,p,l] + dlat * fn_const_dlat * np.cos(lat_rad[l+1]) * 0.25*(np.sum(w_val))	  

            for l in range(len(data['lat'])-1,0,-1 ):  #43-1 includive where there are 44 elems
                w_val_rev = np.array([omega[t,p,l],omega[t,p+1,l-1],omega[t,p,l-1],omega[t,p+1,l]])
              
                dlat=lat_rad[l] - lat_rad[l-1]   #positive
                psi_lat_rev[t,p,l-1] = psi_lat_rev[t,p,l]+ dlat*fn_const_dlat*np.cos(lat_rad[l-1])* 0.25*(np.sum(w_val_rev))
               
	    
    #integrate  bottom to top and then top to bottom [v] 
    print '... integrate second direction'  
    for t in xrange(len(data['time'])): 
        for l in xrange(len(data['lat']) -1):
            fn_const_dp = (2*np.pi*rad_earth*np.cos(lat_rad[l+1]))/grav
            psi_p[t,0,l] = 0
            psi_p_rev[t,0,l] = 0
      
            for p in xrange(len(data['pfull'])-1):
                v_val=np.array([vcomp[t,p,l],vcomp[t,p,l+1],vcomp[t,p+1,l],vcomp[t,p+1,l+1]])
                dp = (data['pfull'][p]-data['pfull'][p+1])*100   #in Pa and negative
                psi_p[t,p+1,l]=psi_p[t,p,l]+dp * fn_const_dp * 0.25*(np.sum(v_val))
	  
            for p in range(len(data['pfull'])-1,0,-1 ):
                v_val_rev=np.array([vcomp[t,p-1,l],vcomp[t,p-1,l+1],vcomp[t,p,l+1],vcomp[t,p,l]])
                dp=(data['pfull'][p] - data['pfull'][p-1])*100   #in Pa and positive
                psi_p_rev[t,p-1,l] = psi_p_rev[t,p,l] + dp*fn_const_dp* 0.25*(np.sum(v_val_rev))
               
	    
    print '  Average the stream functions'   	    
    psi_final = np.zeros([4,len(data['time']),len(data['pfull']),len(data['lat'])])
    psi_final[0,:,:,:] = psi_lat
    psi_final[1,:,:,:] = psi_p
    psi_final[2,:,:,:] = psi_p_rev
    psi_final[3,:,:,:] = psi_lat_rev
  
   #take the mean over the four stream funstions
    psi = MeanOverDim(data=psi_final,dim=0)
  

    return psi

