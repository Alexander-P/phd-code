    
import numpy as np
import matplotlib.pyplot as plt

def plot_twocol(exps, s_slice=None, season=None):

    if season:
        season_array = [np.arange(8,26), np.arange(26,44), np.arange(44,62), np.append(np.arange(0,8),np.arange(62,72))]
        s_sclice = season_array[s_slice]
    tlev = np.arange(180, 320, 20)
    psilev = psilev = np.array([-15e10, -12e10, -9e10, -6e10,-5e10, -4e10, -3e10, -2e10, -1e10,
                                1e10, 2e10, 3e10, 4e10,5e10, 6e10, 9e10, 12e10, 15e10])
    ulev = np.arange(-75, 76, 10)
    
    fig, axes = plt.subplots(len(exps), 2, figsize=(12, 4*len(exps)))        
    for i, exp in enumerate(exps):
        ds = open_climatology(data_dir, exp)
        if season == None:
            dss = ds.mean('pentad')
        else:
            dss = ds.isel(pentad=s_slice).mean('pentad')
            
        x = dss.lat.values
        y = dss.pfull.values    
        X, Y = np.meshgrid(x, y)
        
        axes[i,0].contour(X, Y, dss.ucomp, ulev, add_colorbar=False, robust=True, colors='k')
        dss.psi.plot.contourf(x='lat', y='pfull', ax=axes[i,0], add_colorbar=False, robust=True, levels=psilev)
        axes[i,0].invert_yaxis()
        axes[i,0].set_xlabel('Latitude')
        axes[i,0].set_ylabel('Pressure')
        
        
        dss.temp.plot.contourf(x='lat', y='pfull', ax=axes[i,1], add_colorbar=False, levels=tlev)
        axes[i,1].contour(X, Y, dss.uv, add_colorbar=False, robust=True,  colors='k')
        axes[i,1].invert_yaxis()