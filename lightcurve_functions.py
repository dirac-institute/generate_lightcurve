from __future__ import print_function
import pylab as plt

import string
import os
import numpy as np
import warnings
import re
import argparse
import sys
import glob
import glob
import argparse
import pyslalib.slalib as sla
import random


au_to_meters = 149597870700.0

def compute_oorb_astroidcentric_helio_and_toppo_vectors_with_JD(oorb_location, in_orbit_file_des, start_date_mjd, end_date_mjd, time_step_days, orbit_file_name):
    id_8 = id_generator()
    oorb_command_prop = oorb_location + " --task=propagation --orb-in=" + in_orbit_file_des
    oorb_command_prop_complete = oorb_command_prop + " " + "--epoch-mjd-tt=" + str(start_date_mjd) + " " + "> " + "in_orb" + id_8 + ".des"
    os.system(oorb_command_prop_complete)
    duration = end_date_mjd - start_date_mjd
    oorb_command_ephem = oorb_location + " --task=ephemeris --code=500 --orb-in=in_orb" + id_8 + ".des " + "--timespan=" + str(duration) + " " + "--step=" + str(time_step_days) + " > ephem_" + id_8
    os.system(oorb_command_ephem)
    load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z = np.loadtxt("ephem_" + id_8, usecols=(28,29,30,34,35,36,2,3,9))
    obs_x, obs_y, obs_z = load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z[:,3], load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z[:,4], load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z[:,5]
    hel_x, hel_y, hel_z = load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z[:,0], load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z[:,1], load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z[:,2]
    toppo_x, toppo_y, toppo_z = hel_x - obs_x, hel_y - obs_y, hel_z - obs_z
    astro_toppo_x, astro_toppo_y, astro_toppo_z = -1.0 * toppo_x, -1.0 * toppo_y, -1.0 * toppo_z
    astro_hel_x, astro_hel_y, astro_hel_z = -1.0 * hel_x, -1.0 * hel_y, -1.0 * hel_z
    JD = load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z[:,6] + 2400000.5
    delta_au = load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z[:,7]
    #mags
    m = load_orbit_helio_x_helio_y_helio_z_obs_x_obs_y_obs_z[:,8]
    #light_time_correction
    delta_m = delta_au * au_to_meters
    light_time_travel_days = (delta_m/3e8) / (24*3600.)
    JD_light_time_corrected = JD - light_time_travel_days
    #os.system('rm *' + id_8)
    return np.vstack((JD_light_time_corrected, m, astro_hel_x, astro_hel_y, astro_hel_z, astro_toppo_x, astro_toppo_y, astro_toppo_z)).T

def get_rates(rate, pa, dec_deg, dec_min, dec_sec): #rate is in "/min, pa is in degs
    #print (np.sign(dec_deg) * (dec_min/60.))
    RA = (rate * (1000./60.) * np.sin(np.radians(pa)))/np.cos(np.radians(dec_deg + ((np.sign(dec_deg) * dec_min)/60.) + ((np.sign(dec_deg) * dec_sec)/3600.)))
    DEC = (rate * (1000./60.)) * np.cos(np.radians(pa))
    return RA, DEC #mili arcsec per sec

def get_rates_no_cos_dec(rate, pa_deg): #rate is in "/min, pa is in degs USE FOR UH 88" when cos dec is turned on
    #print (np.sign(dec_deg) * (dec_min/60.))                                                                                                                                                              
    RA = (rate * (1000./60.) * np.sin(np.radians(pa_deg)))
    DEC = (rate * (1000./60.)) * np.cos(np.radians(pa_deg))
    return RA, DEC #mili arcsec per sec  

def get_rates_cos_dec(rate, pa, dec_deg, dec_min, dec_sec): #rate is in "/min, pa is in degs                                                                                                                    
    #print (np.sign(dec_deg) * (dec_min/60.))                                                                                                                                                           
    RA = (rate * (1000./60.) * np.sin(np.radians(pa)))*np.cos(np.radians(dec_deg + ((np.sign(dec_deg) * dec_min)/60.) + ((np.sign(dec_deg) * dec_sec)/3600.)))
    DEC = (rate * (1000./60.)) * np.cos(np.radians(pa))
    return RA, DEC #mili arcsec per sec

def grep_asteroid_from_MPCORBDAT_to_KEP_DES_format(asteroid_numbered_name, mpc_orb_dat_location):

    '''
    grep_asteroid_from_MPCORBDAT_to_KEP_DES_format('808','/Users/bolin/Thermal/asteroid_lists/MPCORB.DAT')
    :param asteroid_numbered_name:
    :param mpc_orb_dat_location: =  /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT
    :return:
    '''

    id_8 = id_generator()
    des_suffix = '1 6 -1 OpenOrb'
    kep_string = 'KEP'
    awk_command = '''awk '{print $'1', $'11', $'9', $'8', $'7', $'6', $'5', $'2', $'24'}' '''
    space = " "

    asteroid_name_paren = '(' + asteroid_numbered_name + ')'
    grep_command_asteroid = '''grep "''' + asteroid_name_paren + '''"''' + " " + mpc_orb_dat_location + " > temp_grep_" + id_8
    os.system(grep_command_asteroid)
    os.system(awk_command + ' < ' + 'temp_grep_'+ id_8 + ' > temp_awk_' + id_8)
    #number_a_au_e_i_deg_Omega_deg_omega_deg_M_deg_H_epoch_YMD
    awk_output = np.loadtxt('temp_awk_' + id_8, dtype='string')
    year = float(awk_output[8][:4])
    month = float(awk_output[8][4:6])
    day = float(awk_output[8][6:])
    mjd = sla.sla_caldj(year, month, day)[0]
    des_out = awk_output[0] + " KEP " + awk_output[1] + space + awk_output[2]  + space + awk_output[3] + space + awk_output[4] + space + awk_output[5] + space + awk_output[6] + space + awk_output[7] + space + str(mjd) + space + des_suffix
    os.system('rm *' + id_8)
    return des_out

def id_generator(size=8, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def M_anom_and_mean_motion_to_time_of_peri(M,n,epoch_mjd_ut): #both M and n are in degrees and degrees/day.
    if M > 180:
        time_peri = epoch_mjd_ut - ((M - 360.) / n)
    if M <=180:
        time_peri = epoch_mjd_ut - (M / n)
    return time_peri

def run_lightcurve_code(JD_helxyz_obs_xyz_array, asteroid_name, shape_model_directory, lightcurve_code, start_time_mjd, end_time_mjd):
    id_generator_lc = id_generator()
    result_file_name = asteroid_name + '_lc_' + str(int(start_time_mjd)) +'_to_' + str(int(end_time_mjd))  + '.txt'

    if len(JD_helxyz_obs_xyz_array)< 1001:
        lc_file_name = 'lc'+ id_generator_lc +'.txt'
        lc_file_name_tmp = 'lc'+ id_generator_lc +'.txttmp'
        np.savetxt(lc_file_name,JD_helxyz_obs_xyz_array)
        os.system('echo 1' + ' > ' + lc_file_name_tmp)
        os.system('echo ' + str(len(JD_helxyz_obs_xyz_array)) + ' 0 >> '+ lc_file_name_tmp)
        os.system('cat ' + lc_file_name + ' >> ' + lc_file_name_tmp)
    if len(JD_helxyz_obs_xyz_array)> 1001: #thousand data point limit per lightcurve
        fraction_left, times_around = np.modf(len(JD_helxyz_obs_xyz_array)/1000.)
        lc_file_name = 'lc'+ id_generator_lc +'.txt'
        lc_file_name_tmp = 'lc'+ id_generator_lc +'.txttmp'
        lc_file_name_array_temp = 'lc'+ id_generator_lc +'.txttmparray'
        os.system('echo ' + str(int(times_around)+1) + ' > ' + lc_file_name_tmp)
        for i in range(0, int(times_around)):
            np.savetxt(lc_file_name_array_temp,JD_helxyz_obs_xyz_array[int((i)*1000):int((i+1)*1000)])
            os.system('echo ' + str(len(JD_helxyz_obs_xyz_array[int((i)*1000):int((i+1)*1000)])) + ' 0 >> '+ lc_file_name_tmp)
            os.system('cat ' + lc_file_name_array_temp + ' >> ' +  lc_file_name_tmp)
        #remainder
        if fraction_left > 0.0:
            np.savetxt(lc_file_name_array_temp,JD_helxyz_obs_xyz_array[int((i+1)*1000):])
            os.system('echo ' + str(len(JD_helxyz_obs_xyz_array[int((i+1)*1000):])) + ' 0 >> '+ lc_file_name_tmp)
            os.system('cat ' + lc_file_name_array_temp + ' >> ' +  lc_file_name_tmp)
    shape_mode_asteroid_directory = shape_model_directory + '/' + asteroid_name + '/'

    lcgenerator_command = 'cat ' + lc_file_name_tmp + ' | ' + lightcurve_code + ' -v ' + shape_mode_asteroid_directory +  '*.spin.txt ' + shape_mode_asteroid_directory + '*.shape.txt ' + result_file_name
    os.system(lcgenerator_command)
    intensities = np.loadtxt(result_file_name)
    (JD_helxyz_obs_xyz_array[:,0], intensities)
    JD_intensity = np.vstack((JD_helxyz_obs_xyz_array[:,0], intensities)).T
    np.savetxt(result_file_name ,JD_intensity)

def string_seperated_to_array_spaces(input_array,data_type_str):
    input_array = ",".join(input_array)
    cols_index_storage = np.array([])
    starting_point = 0
    comma_count = 0
    for i in range(0, len(input_array)):
        test = input_array[starting_point:].find(',')
        #print i, test, starting_point
        if test != -1 and test!=0:
            cols_index_storage = np.append(cols_index_storage,float(input_array[starting_point:input_array[starting_point:].find(',')+starting_point]))
            starting_point = input_array[starting_point:].find(',')+starting_point
        if test == -1:
            cols_index_storage = np.append(cols_index_storage,float(input_array[starting_point:]))
            break
        if test == 0:
            starting_point += 1
            comma_count +=1
    return cols_index_storage.astype(data_type_str)
