#Bryce Bolin
import random
import argparse
import pyslalib.slalib as sla
import sys
sys.path.insert(0, '/Users/bolin/NEO/lightcurves/scripts/')
from lightcurve_functions import *

'''
sample execution:

test

ipython -i -- make_synthetic_lightcurve.py -lcc /Users/bolin/NEO/lightcurves/version_0.2.1/lcgenerator/lcgenerator -smd /Users/bolin/NEO/lightcurves/version_0.2.1/shapemodels -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 306 -stetss 57303 57320 0.0025

ipython -i -- make_synthetic_lightcurve.py -lcc /Users/bolin/NEO/lightcurves/version_0.2.1/lcgenerator/lcgenerator -smd /Users/bolin/NEO/lightcurves/version_0.2.1/shapemodels -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 306 -stetss 44113.0 44114 0.002

ipython -i -- make_synthetic_lightcurve.py -lcc /Users/bolin/NEO/lightcurves/version_0.2.1/lcgenerator/lcgenerator -smd /Users/bolin/NEO/lightcurves/version_0.2.1/shapemodels -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 3200 -stetss 59323 59345 0.0025

ipython -i -- make_synthetic_lightcurve.py -lcc /Users/bolin/NEO/lightcurves/version_0.2.1/lcgenerator/lcgenerator -smd /Users/bolin/NEO/lightcurves/version_0.2.1/shapemodels -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 3200 -stetss 49657.875099 49658.875099 0.0025

must be careful you dont require the generation of more than 10,000 lightcurve datapoints. lcgenerator can't handle more than 10,000.

'''

parser = argparse.ArgumentParser()
parser.add_argument("-lcc", "--lightcurve_code", help="location of lightcurve code, e.g., /Users/bolin/NEO/lightcurves/version_0.2.1/lcgenerator/lcgenerator", nargs='*')
parser.add_argument("-smd", "--shape_model_directory", help="location of lightcurve shape models, e.g., /Users/bolin/NEO/lightcurves/version_0.2.1/shapemodels", nargs='*')
parser.add_argument("-mpcorbf", "--mpc_orb_file", help="location of lightcurve shape models, e.g., /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT", nargs='*')
parser.add_argument("-ooloc", "--oorb_location", help="location of oorb, e.g., /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb", nargs='*')
parser.add_argument("-an", "--asteroid_name", help="numbered name of asteriod", nargs='*')
parser.add_argument("-stetss", "--start_time_end_time_step_size", help="start_time,end_time in MJD and step size in days 57303 57335 0.0025", nargs='*')
args = parser.parse_args()

lightcurve_code = str(args.lightcurve_code[0])
shape_model_directory = str(args.shape_model_directory[0])
mpc_orb_file = str(args.mpc_orb_file[0])
oorb_location = str(args.oorb_location[0])
asteroid_name = str(args.asteroid_name[0])
start_time_mjd, end_time_mjd, step_size_days = string_seperated_to_array_spaces(args.start_time_end_time_step_size,'float')

orbit_line =  grep_asteroid_from_MPCORBDAT_to_KEP_DES_format(asteroid_name,mpc_orb_file)
id_generator_orb = id_generator()
orbit_file_name = '''orbit_''' + id_generator_orb + '''.des'''
echo_orbit_line_to_des = '''echo "''' + orbit_line + '''"''' + ''' > ''' + orbit_file_name
os.system(echo_orbit_line_to_des)
#os.system('rm *' + id_generator_orb + '*')

JD_light_time_corrected_m_astro_hel_x_au_astro_hel_y_au_astro_hel_z_au_astro_toppo_x_au_astro_toppo_y_au_astro_toppo_z = compute_oorb_astroidcentric_helio_and_toppo_vectors_with_JD(oorb_location,orbit_file_name, start_time_mjd, end_time_mjd, step_size_days, orbit_file_name)

run_lightcurve_code(JD_light_time_corrected_m_astro_hel_x_au_astro_hel_y_au_astro_hel_z_au_astro_toppo_x_au_astro_toppo_y_au_astro_toppo_z, asteroid_name, shape_model_directory, lightcurve_code, start_time_mjd, end_time_mjd)
os.system('rm *' + id_generator_orb + '*')

'''
JD_light_time_corrected_m_astro_hel_x_au_astro_hel_y_au_astro_hel_z_au_astro_toppo_x_au_astro_toppo_y_au_astro_toppo_z[:,2:]

import matplotlib.pyplot as plt
plt.ion()
plt.figure()
#JD_intensity306 = np.loadtxt('306_lc_57303_to_57320.txt')
JD_intensity306 = np.loadtxt('306_lc_44113_to_44114.txt')
plt.plot(JD_intensity306[:,0], JD_intensity306[:,1])
plt.xlabel('JD')
plt.ylabel('Relative intensity')
plt.xlim(3.5+2.44411e6,4+2.44411e6)
plt.savefig('306_lc_57309_5_to_57310_0.png')

import matplotlib.pyplot as plt
plt.ion()
plt.figure()
#JD_intensity3200 = np.loadtxt('3200_lc_59323_to_59345.txt')
JD_intensity3200 = np.loadtxt('3200_lc_49657_to_49658.txt')
plt.plot(JD_intensity3200[:,0], JD_intensity3200[:,1])
plt.xlabel('JD')
plt.ylabel('Relative intensity')
plt.xlim(8.5+2.44965e6,9.0+2.44965e6)
plt.savefig('3200_lc_49658_to_49658_5.png')


'''