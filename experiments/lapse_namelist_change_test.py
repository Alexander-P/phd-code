import numpy as np
import os

from gfdl.experiment import Experiment, DiagTable

exp = Experiment('lapse_namelist_old', overwrite_data=True,
	repo='git@github.com:Alexander-P/Isca.git',
    commit='9cb6505262fc4c3b206db6d9ba3c85ef271ba9eb')

exp2 = Experiment('lapse_namelist_new', overwrite_data=True,
	repo='git@github.com:ExeClim/Isca.git',
    commit='bde46187187b4624aea2ea8174a29c2af88c34e0')
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

diag.add_field('dynamics', 'ucomp_sq', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'ucomp_vcomp', time_avg=True, files=['atmos_daily'])
diag.add_field('dynamics', 'ucomp_omega', time_avg=True, files=['atmos_daily'])

diag.add_field('hs_forcing', 'teq', time_avg=True, files=['atmos_daily'])
diag.add_field('hs_forcing', 'h_trop', time_avg=True, files=['atmos_daily'])

exp.use_diag_table(diag)

#Turn off the full, slow radiation scheme compilation
exp.disable_rrtm()

exp.compile()
exp.clear_rundir()

exp.namelist['spectral_dynamics_nml'] = {
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

exp.namelist['main_nml'] = {
    'dt_atmos': 300,
    'days': 90,
    'calendar': 'no_calendar'
}

exp.namelist['atmosphere_nml']['idealized_moist_model'] = False

exp.namelist['hs_forcing_nml'] = {
    'equilibrium_t_option' : 'top_down',
    'ml_depth': 10.,
    'spinup_time': 108000,
    'ka': -20.,
    'ks': -5.
 }

exp.namelist['constants_nml'] = {
    'orbital_period' : 360
 }
 
exp.namelist['astronomy_nml'] = {
    'obliq' : 30
 }
#s End namelist changes from default values


exp2.use_diag_table(diag)
exp2.disable_rrtm()
exp2.compile()
exp2.clear_rundir()
exp2.namelist = exp.namelist.copy()
exp2.namelist['hs_forcing_nml']['lapse_rate'] = 0.666342857
exp2.inputfiles = exp.inputfiles
exp2.execdir = exp.execdir

exp.runmonth(1, use_restart=False, num_cores=8)
exp2.runmonth(1, use_restart=False, num_cores=8)

		


