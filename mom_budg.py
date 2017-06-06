"""
Evaluate and plot momentum budget at 150 hPa

"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#from data_handling import time_means, month_dic
import sh
import gradients as gr
from pylab import rcParams
from useful import *


def partition_advection(data, lons, lev=150):
    
    #First do uu terms
    uu_trans_dx = -86400. * gr.ddx( (data.ucomp_sq - data.ucomp**2).sel(pfull=lev,method='nearest') ) # <u'u'> = <uu> - <u><u>
    
    u = data.ucomp.sel(pfull=lev,method='nearest') # u
    u_dx = -86400. * gr.ddx( u )  # dudx
    
    u_dudx_zav = u.sel(lon=lons,method='nearest').mean('lon') * u_dx.sel(lon=lons,method='nearest').mean('lon') # [u][dudx]

    u_dudx_stat = (u * u_dx).sel(lon=lons,method='nearest').mean('lon') - u_dudx_zav           # u*dudx* = [ududx] - [u][dudx]
    
    data['uu_trans_dx'] = (('xofyear','lat'), uu_trans_dx.sel(lon=lons,method='nearest').mean('lon') )	
    data['u_dudx_stat'] = (('xofyear','lat'), u_dudx_stat )	
    data['u_dudx_zav']  = (('xofyear','lat'), u_dudx_zav )
    
    print 'uu terms done'
    
    #Next do uv terms
    uv_trans_dy = -86400. * gr.ddy( (data.ucomp_vcomp - data.ucomp * data.vcomp).sel(pfull=lev,method='nearest') , uv=True)

    v = data.vcomp.sel(pfull=lev,method='nearest').load() # v
    u_dy = -86400. * gr.ddy( u )  # dudy
    
    v_dudy_zav = v.sel(lon=lons,method='nearest').mean('lon') * u_dy.sel(lon=lons,method='nearest').mean('lon') # [v][dudy]
        
    v_dudy_stat = (v * u_dy).sel(lon=lons,method='nearest').mean('lon') - v_dudy_zav           # v*dudy* = [vdudy] - [v][dudy]
        
    data['uv_trans_dy'] = (('xofyear','lat'), uv_trans_dy.sel(lon=lons,method='nearest').mean('lon'))	
    data['v_dudy_stat'] = (('xofyear','lat'), v_dudy_stat)	
    data['v_dudy_zav']  = (('xofyear','lat'), v_dudy_zav )
    
    print 'uv terms done'
    
    #Finally do uw terms
    uw_trans_dp = -86400. * gr.ddp( (data.ucomp_omega - data.ucomp * data.omega).sel(lon=lons,method='nearest').mean('lon') )
    
    w = data.omega.sel(pfull=lev,method='nearest').load() # w
    u_dp = -86400. * (gr.ddp(data.ucomp)).sel(pfull=lev,method='nearest')  # dudp
    
    w_dudp_zav = w.sel(lon=lons,method='nearest').mean('lon') * u_dp.sel(lon=lons,method='nearest').mean('lon')
    w_dudp_stat = (w * u_dp).sel(lon=lons,method='nearest').mean('lon') - w_dudp_zav
    
    data['uw_trans_dp'] = (('xofyear','lat'), uw_trans_dp.sel(pfull=lev,method='nearest'))	
    data['w_dudp_stat'] = (('xofyear','lat'), w_dudp_stat)	
    data['w_dudp_zav']  = (('xofyear','lat'), w_dudp_zav )	
    
    print 'uw terms done'
    

def prepare_data(data, exp):
    # add required terms to dataset for use with the rest of Ruth's code
    # currently adding ucomp_ucomp, ucomp_vcomp, ucomp_omega
    # returning a climatology
    print('Preparing data')
    
    
    #data['ucomp_ucomp'] = (('time','lat'), data.ucomp*data.ucomp)
    #print('ucomp_ucomp calculated')
    #data['ucomp_vcomp'] = (('time','lat'), data.ucomp*data.vcomp)
    #print('ucomp_vcomp calculated')
    #data['ucomp_omega'] = (('time','lat'), data.ucomp*data.omega)
    #print('ucomp_omega calculated')
    
    data.time.data = data.time.data/24 - 1
    ds = open_zmean_runset('/scratch/ap587/dry_data/', exp)
    ds = ds.isel(time=slice(720,None))
    lag = calc_seasonal_lag(ds.teq)
    data.coords['xofyear'] = np.mod( data.time+45+lag, 360.) //90 + 1
    
    #needs the processed file to exist --- add error handling here
    
    if 'h' in ds.keys():
        data['psi'] = ds.psi
    else:
        data['psi'] = ds.rename({'height':'h'}).psi
    data = data.groupby('xofyear').mean('time')
    
    edges = np.ones([data.xofyear.size,3])
    ### add a term for the edge of the hadley cell locations
    for i in range(data.xofyear.size):
        edges[i,:] = get_hc_edges(data.isel(xofyear=i).mean('lon'), method='psi')

    data['HC_edgeS'] = (('xofyear'), edges[:,0])
    data['HC_edgeE'] = (('xofyear'), edges[:,1])
    data['HC_edgeN'] = (('xofyear'), edges[:,2])
    return data
    
    
    
def mom_budg_hm(exp, lev=150, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 12, 7
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/ap587/plots/momentum/exp2_hc/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = open_runset('/scratch/ap587/dry_data/', exp+'/', range(8,100))
    data = prepare_data(data,exp)
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    #advective terms
    partition_advection(data, lons, lev)
    
    #Coriolis
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    fv = data.vcomp.sel(pfull=lev,method='nearest') * f * 86400.
    fv = fv.sel(lon=lons,method='nearest').mean('lon')
    
    #totp = ((data.convection_rain + data.condensation_rain)*86400.).sel(lon=lons).mean('lon')
    
    abs_vort = (data.vor + f).sel(lon=lons,method='nearest').mean('lon')*86400.
    
    abs_vortv = (abs_vort * data.vcomp.mean('lon')).sel(pfull=lev,method='nearest')    
    vortv = (data.vor * data.vcomp).sel(lon=lons,method='nearest').sel(pfull=lev,method='nearest').mean('lon')*86400.
    
    #Geopotential gradient
    dphidx = gr.ddx(data.height.sel(pfull=lev,method='nearest'))
    dphidx = -86400. * 9.8 * dphidx.sel(lon=lons,method='nearest').mean('lon')
        
    mom_mean = data.u_dudx_zav + data.v_dudy_zav + data.w_dudp_zav
    mom_trans = data.uu_trans_dx + data.uv_trans_dy + data.uw_trans_dp
    mom_stat = data.u_dudx_stat + data.v_dudy_stat + data.w_dudp_stat
    
    mom_sum = fv + dphidx + mom_mean + mom_trans + mom_stat
    vort_eddy_sum = abs_vortv + mom_trans + mom_stat
    
    levels = np.arange(-6,6.1,1.)
    levels_max = max((abs(abs_vortv)).max().values, abs(mom_trans).max().values)
    levels = np.linspace(-np.floor(levels_max), np.floor(levels_max), 12)
    
    #mn_dic = month_dic(1)
    #tickspace = range(13,72,18)
    tickspace = data.xofyear.values
    #labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    labels = ['%.d' % k for k in tickspace]
    
    # Six subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f1=fv.plot.contourf(ax=ax1, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax1, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax1.set_ylabel('Latitude')
    ax1.set_title('Coriolis', fontsize=17)
    ax1.set_ylim(-90,90)
    ax1.grid(True,linestyle=':')
    ax1.set_yticks(np.arange(-60.,61.,30.))
    #ax1.text(-15, 60, 'a)')
    
    #Second plot
    mom_mean.plot.contourf(ax=ax2, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #ax2.contour(data.xofyear, data.lat, abs_vort.sel(pfull=lev).T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2)
    #abs_vort.sel(pfull=lev, method='nearest').plot.contour(ax=ax2, x='xofyear', y='lat', extend='both', levels=levels, add_colorbar=False, add_labels=False, colors='k', linewidths=2)
    #totp.plot.contour(ax=ax2, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax2.grid(True,linestyle=':')
    ax2.set_title('Mean state advection', fontsize=17)
    ax2.set_ylim(-90,90)
    #ax2.text(-5, 60, 'b)')
    
    #Third plot
    dphidx.plot.contourf(ax=ax3, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax3, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax3.grid(True,linestyle=':')
    ax3.set_ylim(-90,90)
    ax3.set_title('Geopotential gradient', fontsize=17)
    #ax3.text(-5, 60, 'c)')
    
    #Fourth plot
    mom_trans.plot.contourf(ax=ax4, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax4, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax4.grid(True,linestyle=':')
    ax4.set_ylabel('Latitude')
    ax4.set_xticks(tickspace)
    ax4.set_xticklabels(labels,rotation=25)
    ax4.set_ylim(-90,90)
    ax4.set_title('Transient eddy flux conv.', fontsize=17)
    ax4.set_yticks(np.arange(-60.,61.,30.))
    #ax4.text(-15, 60, 'd)')
    
    #Fifth plot
    mom_stat.plot.contourf(ax=ax5, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax5, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax5.grid(True,linestyle=':')
    ax5.set_xticks(tickspace)
    ax5.set_xticklabels(labels,rotation=25)
    ax5.set_title('Stat. eddy flux conv.', fontsize=17)
    ax5.set_ylim(-90,90)
    #ax5.text(-5, 60, 'e)')
    
    #Sixth plot
    mom_sum.plot.contourf(ax=ax6, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax6, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax6.grid(True,linestyle=':')
    ax6.set_xticks(tickspace)
    ax6.set_xticklabels(labels,rotation=25)
    ax6.set_ylim(-90,90)
    ax6.set_title('Residual', fontsize=17)
    #ax6.text(-5, 60, 'f)')
    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.25, wspace=0.12)
    #Colorbar
    cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    #cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True,fraction=0.15, aspect=30)
    #cb1.set_label('$ms^{-1}day^{-1}$')
    
    if lonin == [-1.,361.]:
        figname = exp+'_zon_mom_budg_' + str(lev) + 'hPa' + '.pdf'
    else:
        figname = exp+'_zon_mom_budg_' + str(lev) + 'hPa' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    ##########################
    ##########################
    
    # Six subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f1=vort_eddy_sum.plot.contourf(ax=ax1, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax1, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax1.set_ylabel('Latitude')
    ax1.set_title('Residual abs vorticity + transient eddy flux conv.', fontsize=17)
    ax1.set_ylim(-90,90)
    ax1.grid(True,linestyle=':')
    ax1.set_yticks(np.arange(-60.,61.,30.))
    #ax1.text(-15, 60, 'a)')
    
    #Second plot
    fv.plot.contourf(ax=ax2, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #ax2.contour(data.xofyear, data.lat, abs_vort.sel(pfull=lev).T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2)
    #abs_vort.sel(pfull=lev,method='nearest').plot.contour(ax=ax2, x='xofyear', y='lat', extend='both', levels=levels, add_colorbar=False, add_labels=False, colors='k', linewidths=2)
    #totp.plot.contour(ax=ax2, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax2.grid(True,linestyle=':')
    ax2.set_title('Coriolis', fontsize=17)
    ax2.set_ylim(-90,90)
    #ax2.text(-5, 60, 'b)')
    
    #Third plot
    vortv.plot.contourf(ax=ax3, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax3, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax3.grid(True,linestyle=':')
    ax3.set_ylim(-90,90)
    ax3.set_title('Relative vorticity', fontsize=17)
    #ax3.text(-5, 60, 'c)')
    
    #Fourth plot
    data.uu_trans_dx.plot.contourf(ax=ax4, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax4, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax4.grid(True,linestyle=':')
    ax4.set_ylabel('Latitude')
    ax4.set_xticks(tickspace)
    ax4.set_xticklabels(labels,rotation=25)
    ax4.set_ylim(-90,90)
    ax4.set_title('uu transient eddy flux conv.', fontsize=17)
    ax4.set_yticks(np.arange(-60.,61.,30.))
    #ax4.text(-15, 60, 'd)')
    
    #Fifth plot
    data.uv_trans_dy.plot.contourf(ax=ax5, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax5, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax5.grid(True,linestyle=':')
    ax5.set_xticks(tickspace)
    ax5.set_xticklabels(labels,rotation=25)
    ax5.set_title('uv trans. eddy flux conv.', fontsize=17)
    ax5.set_ylim(-90,90)
    #ax5.text(-5, 60, 'e)')
    
    #Sixth plot
    data.uw_trans_dp.plot.contourf(ax=ax6, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #totp.plot.contour(ax=ax6, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax6.grid(True,linestyle=':')
    ax6.set_xticks(tickspace)
    ax6.set_xticklabels(labels,rotation=25)
    ax6.set_ylim(-90,90)
    ax6.set_title('uw trans. eddy flux conv.', fontsize=17)
    #ax6.text(-5, 60, 'f)')
    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.25, wspace=0.12)
    #Colorbar
    cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    #cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True,fraction=0.15, aspect=30)
    #cb1.set_label('$ms^{-1}day^{-1}$')
    
    if lonin == [-1.,361.]:
        figname = exp+'_simple_zon_mom_budg_' + str(lev) + 'hPa' + '.pdf'
    else:
        figname = exp+'_simple_zon_mom_budg_' + str(lev) + 'hPa' + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    ############
    #####    
    
    fig, ax = plt.subplots()
    f1 = vort_eddy_sum.mean('xofyear').plot(ax=ax,label='residual')
    fv.mean('xofyear').plot(ax=ax,label='Coriolis (fv)')
    vortv.mean('xofyear').plot(ax=ax,label='Relative vort (zeta*v)')
    data.uu_trans_dx.mean('xofyear').plot(ax=ax,label="- du'u'/dx")
    data.uv_trans_dy.mean('xofyear').plot(ax=ax,label="- du'v'/dy")
    data.uw_trans_dp.mean('xofyear').plot(ax=ax,label="- du'w'/dp")
    (data.u_dudx_stat + data.v_dudy_stat + data.w_dudp_stat).mean('xofyear').plot(ax=ax,label='stat. eddy')
    plt.axvline(x=data.HC_edgeN.mean('xofyear'), ls='--')
    plt.axvline(x=data.HC_edgeE.mean('xofyear'), ls='--')
    plt.axvline(x=data.HC_edgeS.mean('xofyear'), ls='--')    
    
    ax.set_ylabel('Momentum (m/s/day)')
    ax.set_xlabel('Latitude')
    ax.set_xticks(np.arange(-60,61,30))
    ax.text(data.HC_edgeN.mean('xofyear')+3, min(fv.mean('xofyear'))-0.3, 'Hadley cell edges', color='blue',fontsize=12)
    plt.legend(prop={'size':10})
    
    if lonin == [-1.,361.]:
        figname = exp+'_line_zon_mom_budg_' + str(lev) + 'hPa' + '.pdf'
    else:
        figname = exp+'_line_zon_mom_budg_' + str(lev) + 'hPa' + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    ############
    #####  
    
    for i in range(data.xofyear.size):
        
        fig, ax1 = plt.subplots()
        ax1.yaxis.tick_right()
        levels = [-15e10, -12e10, -9e10, -6e10,-5e10, -4e10, -3e10, -2e10, -1e10, 1e10, 2e10, 3e10, 4e10,5e10, 6e10, 9e10, 12e10, 15e10]
        
        x = data.lat.values
        y = data.pfull.values #np.linspace(-8,8,len(data.pfull.values))
        X,Y = np.meshgrid(x,y)
        psi = data.psi.isel(xofyear=i)
        #psi = np.tile(psi,(len(y),1)) only use if want to fill with one pressure level
        cs = ax1.contourf(X, Y, psi, levels, cmap='RdBu_r', alpha=1, add_colorbar=False, extend='both')
        cs.set_clim([-15e10,15e10])
        ax1.set_ylabel('Pressure (mb)')
        ax1.yaxis.set_label_position('right')
        ax1.invert_yaxis()
        
        ax2 = fig.add_subplot(111,sharex=ax1,frameon=False)
        ax2.yaxis.tick_left()
        f1 = mom_sum.isel(xofyear=i).plot(ax=ax2,color='k',label='residual',lw=3)
        #fv.isel(xofyear=i).plot(ax=ax,label='Coriolis (fv)')
        #vortv.isel(xofyear=i).plot(ax=ax,label='Relative vort (zeta*v)')
        #data.uv_trans_dy.isel(xofyear=i).plot(ax=ax,label="- du'v'/dy")
        #data.uw_trans_dp.isel(xofyear=i).plot(ax=ax,label="- du'w'/dp")
        
        (fv+vortv).isel(xofyear=i).plot(ax=ax2, label='Absolute vorticity', ls='--',color='k',lw=3)
        
        fv.isel(xofyear=i).plot(ax=ax, label='Coriolis', ls='--',color='g',lw=2)
        vortv.isel(xofyear=i).plot(ax=ax, label='Relative vorticity', ls='--',color='r',lw=2)
        
        mom_trans.isel(xofyear=i).plot(ax=ax2, label='Transient eddies',ls=':',color='k', lw=3)
        
        data.uv_trans_dy.isel(xofyear=i).plot(ax=ax,label="- du'v'/dy", ls=':',color='g',lw=2)
        data.uw_trans_dp.isel(xofyear=i).plot(ax=ax,label="- du'w'/dp", ls=':',color='r',lw=2)       
        
        
        #plt.axvline(x=data.HC_edgeN.isel(xofyear=i), ls='--')
        #plt.axvline(x=data.HC_edgeE.isel(xofyear=i), ls='--')
        #plt.axvline(x=data.HC_edgeS.isel(xofyear=i), ls='--')  
        
        ax2.set_title('')
        ax2.set_title(exp + 'season ' + str(i))
        ax2.set_ylabel('')
        ax2.set_ylabel('Momentum change (m/s/day)')
        ax2.set_xlabel('')
        ax2.set_xlabel('Latitude')
        ax2.set_xticks(np.arange(-60,61,30))
        ax2.set_xlim([-90,90])
        ax2.set_ylim([-8,8])
        
        ax2.text(-80,6, '(b)', color='k', fontsize=20)
        #ax.text(data.HC_edgeE.isel(xofyear=i)+3, min(fv.isel(xofyear=i))-0.3, 'Hadley cell edges', color='blue',fontsize=12)
        
        

        plt.legend(prop={'size':10})
    
        if lonin == [-1.,361.]:
            figname = 'poster1_'+exp+'season'+str(i)+'_line_zon_mom_budg_' + str(lev) + 'hPa'
        else:
            figname = 'poster1_'+exp+'season'+str(i)+'_line_zon_mom_budg_' + str(lev) + 'hPa' + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1]))
        try:
            plt.savefig(plot_dir + figname+ '.png', dpi=200, transparent=True,bbox_inches='tight')
        except:
            plt.savefig(plot_dir+figname+'.pdf',dpi=200, transparent=True,bbox_inches='tight')
        plt.savefig(plot_dir + figname + '.pdf', dpi=200, transparent=True,bbox_inches='tight')
        plt.close()
    
    

#mom_budg_hm('exp8-ml10000_tau1.0',lev=200)
exp = 'exp2'
exps = []
vals = []
for name in os.listdir('/scratch/ap587/dry_data/'):
    if name.split('_')[0] == exp:        
        exps += [name]
        vals += [name.split('_')[1]]

#exps = ['exp2_hc50.0', 'exp4_HS_hc2000.0', 'exp5_obl60.0', 'exp8-ml50_tau0.1']
#exps = ['exp4_HS_hc2000.0']
for exp in exps:
    print(exp)
    mom_budg_hm(exp,lev=200)