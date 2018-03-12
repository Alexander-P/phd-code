import numpy as np
import os
from gfdl.experiment import Experiment, DiagTable
import f90nml

#Define our base experiment to compile
base_dir=os.getcwd()
GFDL_BASE        = os.environ['GFDL_BASE']

baseexp = Experiment('rrtm', overwrite_data=False)

baseexp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_daily', 24, 'hours')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'ucomp_sq', time_avg=True)
diag.add_field('dynamics', 'ucomp_vcomp', time_avg=True)
diag.add_field('dynamics', 'ucomp_omega', time_avg=True)

baseexp.use_diag_table(diag)


#Compile model if not already compiled
baseexp.compile()

#Empty the run directory ready to run
baseexp.clear_rundir()

#Define values for the 'core' namelist
baseexp.namelist['main_nml'] = f90nml.Namelist({
    'days'   : 30,
    'hours'  : 0,
    'minutes': 0,
    'seconds': 0,
    'dt_atmos':720,
    'current_date' : [1,1,1,0,0,0],
    'calendar' : 'thirty_day'
})


baseexp.namelist['mixed_layer_nml']['depth'] = 10                          #Depth of mixed layer used
baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.31                  #Albedo value used

baseexp.namelist['spectral_dynamics_nml']['graceful_shutdown'] = False
baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 25               #How many model pressure levels to use
baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = False     #Use grey radiation
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = True
baseexp.namelist['idealized_moist_phys_nml']['mixed_layer_bc'] = True
baseexp.namelist['idealized_moist_phys_nml']['convection_scheme'] = 'SIMPLE_BETTS_MILLER' #Use the simple Betts Miller convection scheme from Frierson

baseexp.namelist['spectral_dynamics_nml']['vert_coord_option'] = 'input'   #Use the vertical levels from Frierson 2006
baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 5000.           #Bottom of the model's sponge down to 50hPa (units are Pa)
baseexp.namelist['damping_driver_nml']['trayfric'] = -0.25                 #Drag timescale for model's sponge

# Turn on vertical diffusion in the free atmosphere (recommended for axisymmetric)
baseexp.namelist['vert_turb_driver_nml']['do_diffusivity'] = True         #
baseexp.namelist['vert_turb_driver_nml']['do_mellor_yamada'] = False    #
baseexp.namelist['vert_turb_driver_nml']['use_tau'] =False              # these are set this way
baseexp.namelist['vert_turb_driver_nml']['constant_gust'] = 0.0         # by default in phys.nml
baseexp.namelist['idealized_moist_phys_nml']['turb'] = True				#
##
baseexp.namelist['diffusivity_nml']['free_atm_diff'] = False
##

baseexp.namelist['rrtm_radiation_nml']['solr_cnst'] = 1360.
baseexp.namelist['rrtm_radiation_nml']['dt_rad'] = 4320

obls = np.arange(0, 91, 15)

for obl in obls:
	exp = baseexp.derive('rrtm-obl{}'.format(obl))
	exp.clear_rundir()
	exp.namelist['astronomy_nml']['obliq'] = obl
	exp.screen_runmonth_prefix = 'rrtm obl %.1f' % (obl)
	exp.runmonth(1, use_restart=False,num_cores=16, light=False)
	for i in range(2,33):
		exp.runmonth(i, num_cores=16, light=False)
    	

