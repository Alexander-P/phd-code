import numpy as np
import os
from gfdl.experiment import Experiment, DiagTable
import f90nml

#Define our base experiment to compile
base_dir=os.getcwd()
GFDL_BASE        = os.environ['GFDL_BASE']

baseexp = Experiment('frierson-axi', overwrite_data=False)

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

#Turn off the full, slow radiation scheme compilation

baseexp.disable_rrtm()

#Compile model if not already compiled
baseexp.compile()

#Empty the run directory ready to run
baseexp.clear_rundir()

#Define values for the 'core' namelist
baseexp.namelist['main_nml'] = f90nml.Namelist({
     'days'   : 90,
     'dt_atmos':720,
     'calendar': 'no_calendar'
})


baseexp.namelist['mixed_layer_nml']['depth'] = 10                          #Depth of mixed layer used
baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.31                  #Albedo value used


baseexp.namelist['spectral_dynamics_nml']['num_fourier'] = 2
baseexp.namelist['spectral_dynamics_nml']['lon_max'] = 8
baseexp.namelist['spectral_dynamics_nml']['make_symmetric'] = True
baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 25
baseexp.namelist['spectral_dynamics_nml']['vert_coord_option'] = 'input'

baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True     #Use grey radiation
baseexp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'frierson'        #Select radiation scheme to use, which in this case is Frierson
baseexp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True          #do_seasonal=false uses the p2 insolation profile from Frierson 2006. do_seasonal=True uses the GFDL astronomy module to calculate seasonally-varying insolation.
baseexp.namelist['two_stream_gray_rad_nml']['use_time_average_coszen'] = False
baseexp.namelist['idealized_moist_phys_nml']['convection_scheme'] ='SIMPLE_BETTS_MILLER' #Use the simple Betts Miller convection scheme from Frierson

baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 5000.           #Bottom of the model's sponge down to 50hPa (units are Pa)
baseexp.namelist['damping_driver_nml']['trayfric'] = -0.25                 #Drag timescale for model's sponge

# Turn on vertical diffusion in the free atmosphere (recommended for axisymmetric)
baseexp.namelist['vert_turb_driver_nml']['do_diffusivity'] = True         #
baseexp.namelist['vert_turb_driver_nml']['do_mellor_yamada'] = False    #
baseexp.namelist['vert_turb_driver_nml']['use_tau'] =False              # these are set this way
baseexp.namelist['vert_turb_driver_nml']['constant_gust'] = 0.0         # by default in phys.nml
baseexp.namelist['idealized_moist_phys_nml']['turb'] = True				#
##
baseexp.namelist['diffusivity_nml']['free_atm_diff'] = True
##

obls = np.arange(0, 91, 15)

for obl in obls:
	exp = baseexp.derive('frierson-axi-diffusion_obl{}'.format(obl))
	exp.clear_rundir()
	exp.namelist['astronomy_nml']['obliq'] = obl
	exp.screen_runmonth_prefix = 'frierson obl %.1f' % (obl)
	exp.runmonth(1, use_restart=False,num_cores=2, light=False)
	for i in range(2,33):
		exp.runmonth(i, num_cores=2, light=False)


