#Bryce Bolin
import sys
sys.path.insert(0, '/Users/bolin/NEO/lightcurves/scripts/')
from lightcurve_functions import *

'''
sample execution:

TO DO :

loopify the >10,000 pts

generate 306 2 months 30 second intervals

test
#306
ipython -i -- make_synthetic_lightcurve.py -an 306 -stetss 44113.0 44114 0.002

#487
ipython -i -- make_synthetic_lightcurve.py -an 487 -stetss 44113.0 44124.0 0.002

#890
ipython -i -- make_synthetic_lightcurve.py -an 890 -stetss 44113.0 44124.0 0.002

#1388
ipython -i -- make_synthetic_lightcurve.py -an 1388 -stetss 46313.0 46324.0 0.002

#1723
ipython -i -- make_synthetic_lightcurve.py -an 1723 -stetss 46313.0 46324.0 0.002

#4077
ipython -i -- make_synthetic_lightcurve.py -an 4077 -stetss 44113.0 44120.0 0.001

#4800

ipython -i -- make_synthetic_lightcurve.py -an 4800 -stetss 44113.0 44120.0 0.001

#6136

ipython -i -- make_synthetic_lightcurve.py -an 6136 -stetss 44113.0 44120.0 0.001

#for paper

#1291
ipython -i -- make_synthetic_lightcurve.py -an 1291 -stetss 44113.0 44120.0 0.001

#1388
ipython -i -- make_synthetic_lightcurve.py -an 1388 -stetss 44113.0 44120.0 0.001

#3200
ipython -i -- make_synthetic_lightcurve.py -an 3200 -stetss 49657 49668 0.0025

must be careful you dont require the generation of more than 10,000 lightcurve datapoints. lcgenerator can't handle more than 10,000.

'''

parser = argparse.ArgumentParser()
parser.add_argument("-an", "--asteroid_name", help="numbered name of asteriod", nargs='*')
parser.add_argument("-stetss", "--start_time_end_time_step_size", help="start_time,end_time in MJD and step size in days 57303 57335 0.0025", nargs='*')
args = parser.parse_args()
asteroid_name = str(args.asteroid_name[0])
start_time_mjd, end_time_mjd, step_size_days = string_seperated_to_array_spaces(args.start_time_end_time_step_size,'float')

#set up directories and located necessary files
if not os.path.isdir('./shapemodels/'):os.system('mkdir shapemodels')
shape_model_directory='./shapemodels/'
if not os.path.exists("./shapemodels/lcgenloc"):os.system('mdfind -name lcgenerator | grep "lcgenerator/lcgenerator" | head -1 > shapemodels/lcgenloc')
lightcurve_code = str(np.loadtxt('shapemodels/lcgenloc',dtype='string'))
if not os.path.exists("./shapemodels/mpcorbloc"):os.system('mdfind -name MPCORB.DAT > shapemodels/mpcorbloc')
mpc_orb_file = str(np.loadtxt('shapemodels/mpcorbloc',dtype='string'))
if not os.path.exists("./shapemodels/oorbloc"):os.system('mdfind -name oorb | grep "oorb-master/main/oorb" | grep -v "f90" > shapemodels/oorbloc')
oorb_location = str(np.loadtxt('shapemodels/oorbloc',dtype='string'))[:str(np.loadtxt('shapemodels/oorbloc',dtype='string')).find('oorb.')+4]
if not os.path.exists("./shapemodels/"):os.system('mdfind -name oorb | grep "oorb-master/main/oorb" | grep -v "f90" > shapemodels/oorbloc')

#download shape models
if not os.path.exists("./shapemodels/db_export_simple.php"): os.system('wget http://astro.troja.mff.cuni.cz/projects/asteroids3D/php/db_export_simple.php; mv db_export_simple.php shapemodels')
model_number_asteroid_name_shape_model_version = np.loadtxt('shapemodels/db_export_simple.php', usecols=(0,1,2), dtype='string')#use perferred version in 2nd column
damit_asteroid_number = model_number_asteroid_name_shape_model_version[np.where(model_number_asteroid_name_shape_model_version[:,1]==asteroid_name)][0][0]
model_number = model_number_asteroid_name_shape_model_version[np.where(model_number_asteroid_name_shape_model_version[:,1]==asteroid_name)][0][-1]
identifier = model_number_asteroid_name_shape_model_version[np.where(model_number_asteroid_name_shape_model_version[:,1]==asteroid_name)][0][0]

if not os.path.exists("shapemodels/db_export_simple.php"):os.system('wget -r -np -nd http://astro.troja.mff.cuni.cz/projects/asteroids3D/php/spin.txt.php?model_id='++'; mv db_export_simple.php shapemodels')


if not os.path.isdir('shapemodels/'+ asteroid_name):
    os.system('mkdir shapemodels/'+asteroid_name)
    os.system('wget -r -np -nd http://astro.troja.mff.cuni.cz/projects/asteroids3D/php/spin.txt.php?model_id=' + model_number + '; mv spin.txt.php?model_id=' + model_number + ' shapemodels/'+asteroid_name)
    if int(damit_asteroid_number) < 1000:
        os.system('wget -r -np -nd http://astro.troja.mff.cuni.cz/projects/asteroids3D/data/archive/1-1000/A'+ damit_asteroid_number+ '.M' + model_number + '.shape.txt; mv A'+ damit_asteroid_number+ '.M' + model_number + '.shape.txt shapemodels/'+ asteroid_name)
    if int(damit_asteroid_number) > 1000:
        os.system('wget -r -np -nd http://astro.troja.mff.cuni.cz/projects/asteroids3D/data/archive/1001-2000/A'+ damit_asteroid_number + '.M' + model_number + '.shape.txt; mv A'+ damit_asteroid_number+ '.M' + model_number + '.shape.txt shapemodels/'+ asteroid_name)

orbit_line =  grep_asteroid_from_MPCORBDAT_to_KEP_DES_format(asteroid_name,str(mpc_orb_file))
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

#306
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity306 = np.loadtxt('306_lc_44113_to_44114.txt')
plt.plot(JD_intensity306[:,0], JD_intensity306[:,1],'.')
plt.xlabel('JD')
plt.ylabel('Relative intensity')
#plt.xlim(3.5+2.44411e6,4+2.44411e6)
#plt.savefig('306_lc_57309_5_to_57310_0.png')

#487
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity306 = np.loadtxt('487_lc_44113_to_44124.txt')
plt.plot(JD_intensity306[:,0], JD_intensity306[:,1],'.')
plt.xlabel('JD')
plt.ylabel('Relative intensity')

#890
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity890 = np.loadtxt('890_lc_44113_to_44124.txt')
plt.plot(JD_intensity890[:,0], JD_intensity890[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')
#plt.savefig('890_lc_44113_to_44124.png')

#1291 2nd one
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity890 = np.loadtxt('1291_lc_44113_to_44120.txt')
plt.plot(JD_intensity890[:,0], JD_intensity890[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')


#1339
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity890 = np.loadtxt('1339_lc_44113_to_44120.txt')
plt.plot(JD_intensity890[:,0], JD_intensity890[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')
#plt.savefig('890_lc_44113_to_44124.png')

#1353
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity890 = np.loadtxt('1353_lc_44113_to_44120.txt')
plt.plot(JD_intensity890[:,0], JD_intensity890[:,1],'.')
plt.xlabel('JD')
plt.ylabel('Relative intensity')


#1364
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity890 = np.loadtxt('1364_lc_44113_to_44120.txt')
plt.plot(JD_intensity890[:,0], JD_intensity890[:,1],'.')
plt.xlabel('JD')
plt.ylabel('Relative intensity')

#1388 best so far
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity890 = np.loadtxt('1388_lc_44113_to_44120.txt')
plt.plot(JD_intensity890[:,0], JD_intensity890[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')

#1557
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity890 = np.loadtxt('1557_lc_44113_to_44120.txt')
plt.plot(JD_intensity890[:,0], JD_intensity890[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')

#1641
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity890 = np.loadtxt('1641_lc_44113_to_44120.txt')
plt.plot(JD_intensity890[:,0], JD_intensity890[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')


#1723
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity1723 = np.loadtxt('1723_lc_46313_to_46324.txt')
plt.plot(JD_intensity1723[:,0], JD_intensity1723[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')
#plt.savefig('1723_lc_46313_to_46324.png')

#4077
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity1723 = np.loadtxt('4077_lc_44113_to_44120.txt')
plt.plot(JD_intensity1723[:,0], JD_intensity1723[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')

#4800
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity1723 = np.loadtxt('4800_lc_44113_to_44120.txt')
plt.plot(JD_intensity1723[:,0], JD_intensity1723[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')

#6136 candidate
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity1723 = np.loadtxt('6136_lc_44113_to_44120.txt')
plt.plot(JD_intensity1723[:,0], JD_intensity1723[:,1],'-')
plt.xlabel('JD')
plt.ylabel('Relative intensity')

#3200
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
JD_intensity3200 = np.loadtxt('3200_lc_49657_to_49668.txt')
plt.plot(JD_intensity3200[:,0], JD_intensity3200[:,1])
plt.xlabel('JD')
plt.ylabel('Relative intensity')
#plt.xlim(8.5+2.44965e6,9.0+2.44965e6)
#plt.savefig('3200_lc_49658_to_49658_5.png')




'''