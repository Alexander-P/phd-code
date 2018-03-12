import numpy as np
import os

from gfdl.experiment import Experiment, DiagTable

baseexp = Experiment('top_down_yearlength', overwrite_data=True,
	repo='git@github.com:ExeClim/Isca.git')
    #commit='765df261f6ecbe7abc4fc3598a5770f7bb02d9e2')

diag = DiagTable()

diag.add_file('atmos_daily', 1, 'days', time_units='days')

# Define diag table entries 

diag.add_field('dynamics', 'ps', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'bk', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'pk', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'vor', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'div', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'ucomp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'vcomp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'temp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'slp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'omega', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'height', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'height_half', time_avg=True, files=['atmos_daily'])

diag.add_field('hs_forcing', 'teq', time_avg=True, files=['atmos_daily'])
diag.add_field('hs_forcing', 'h_trop', time_avg=True, files=['atmos_daily'])

baseexp.use_diag_table(diag)

#Turn off the full, slow radiation scheme compilation

baseexp.disable_rrtm()

baseexp.compile()

baseexp.clear_rundir()

baseexp.namelist['spectral_dynamics_nml'] = {
	'num_levels': 25,
	'exponent': 2.5,
	'scale_heights': 4,
	'surf_res': 0.1,
	'robert_coeff': 4e-2,
	'do_water_correction': False,
	'vert_coord_option': 'even_sigma',
	'initial_sphum': 0.,
	'valid_range_T': [0, 700]
}

baseexp.namelist['main_nml'] = {
    'dt_atmos': 300,
    'days': 90,
    'calendar': 'no_calendar'
}

baseexp.namelist['atmosphere_nml']['idealized_moist_model'] = False


baseexp.namelist['hs_forcing_nml'] = {
    'equilibrium_t_option' : 'top_down',
    'ml_depth': 10.,
    'spinup_time': 108000,
    'ka': -20.,
    'ks': -5.
 }

baseexp.namelist['constants_nml'] = {
    'orbital_period' : 360
 }
 
baseexp.namelist['astronomy_nml'] = {
    'obliq' : 30
 }

yls = [20]
obls = [60]

#s End namelist changes from default values


for yl in yls:
    for obl in obls:
        exp = baseexp.derive('topdown_hc10_yl{}-obl{}'.format(yl,obl))
        exp.clear_rundir()

        exp.use_diag_table(diag)
        exp.execdir = baseexp.execdir

        exp.inputfiles = baseexp.inputfiles

        exp.namelist = baseexp.namelist.copy()
        exp.namelist['astronomy_nml']['obliq'] = obl
        exp.namelist['constants_nml']['orbital_period'] = yl*360
        exp.namelist['hs_forcing_nml']['ml_depth'] = 10
    
        exp.runmonth(1, use_restart=False, num_cores=16)
        for i in range(2, 9*yl):  # 81, 161
            exp.runmonth(i, num_cores=16)
		


