#!/usr/bin/python

import math
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from matplotlib import cm
import matplotlib.colorbar as cb
from scipy.interpolate import interp1d

# SANS merging procedure

def read_plot_NIMROD(path, file):

    pathfile = path+'/'+file
    file_to_read = open(pathfile, 'r')

    nlines = 0

    for line in file_to_read:

        nlines = nlines+1

    file_to_read.seek(0)

    for iosef in range(14):
        osef = file_to_read.readline()

    q_file = np.zeros(nlines-15)
    dcs_file = np.zeros((nlines-15))
    err_file = np.zeros((nlines - 15))

    for iline in range(nlines-15):

        line_to_read = file_to_read.readline()
        line_cut = line_to_read.split()
        q_file[iline] = float(line_cut[0])
        dcs_file[iline] = 0.094*float(line_cut[1])  #0.094*
        err_file[iline] = 0.094*float(line_cut[2])

    file_to_read.close()

    return nlines-15, q_file, dcs_file, err_file

def read_plot_SANS2D(path, file):

    pathfile = path+'/'+file
    file_to_read = open(pathfile, 'r')

    nlines = 0

    for line in file_to_read:

        nlines = nlines+1

    file_to_read.seek(0)
    
    for iosef in range(5):
        osef = file_to_read.readline()


    q_file = np.zeros(nlines-5)
    dcs_file = np.zeros((nlines-5))
    err_file = np.zeros((nlines-5))

    for iline in range(nlines-5):

        line_to_read = file_to_read.readline()
        line_cut = line_to_read.split()
        q_file[iline] = float(line_cut[0])
        dcs_file[iline] = float(line_cut[1]) 
        err_file[iline] = float(line_cut[2]) 

    file_to_read.close()

    return nlines -5, q_file, dcs_file, err_file

path_SANS2D = '../Sans2d'
path_NIMROD = '../NIMROD'

list_path_SANS2D = listdir(path_SANS2D)
list_path_NIMROD = listdir(path_NIMROD)



max_q_SANS2D = 0.4
count = 0

#reading of SANS2D and cut of data when it starts being negative

file_SANS2D = 'FILE.txt'

nlines_SANS2D, q_SANS2D, dcs_SANS2D, err_SANS2D = read_plot_SANS2D(path_SANS2D, file_SANS2D)

for ilines in range(nlines_SANS2D):

    dcs_SANS2D[ilines] = dcs_SANS2D[ilines]

    if dcs_SANS2D[ilines] < 0 and count == 0:

        max_q_SANS2D = q_SANS2D[ilines-1]
        count = 1

    file_NIMROD = 'FILE.mint01'

    nlines_NIMROD, q_NIMROD, dcs_NIMROD, err_NIMROD = read_plot_NIMROD(path_NIMROD, file_NIMROD)

    for ilines in range(nlines_NIMROD):

        dcs_NIMROD[ilines] = dcs_NIMROD[ilines]
        err_NIMROD[ilines] = err_NIMROD[ilines]

# interpolation of NIMROD and SANS2D

f_NIMROD = interp1d(q_NIMROD, dcs_NIMROD)
f_SANS2D_1 = interp1d(q_SANS2D, dcs_SANS2D)

f_NIMROD_err = interp1d(q_NIMROD, err_NIMROD)
f_SANS2D_1_err = interp1d(q_SANS2D, err_SANS2D)

chi_square_1 = np.zeros(2000)

#Construction of the new Q range combining the Q values from NIMROD and SANS2d

q_tot_1 = np.zeros(nlines_NIMROD+nlines_SANS2D)

for i in range(nlines_SANS2D):

    q_tot_1[i] = q_SANS2D[i]

for i in range(nlines_NIMROD):

    q_tot_1[i+nlines_SANS2D] = q_NIMROD[i]

# LEVELLING ; multiplicatio of SANS2D data by an offset and calculation of a chi-square difference with NIMROD's data
# then find optimal offset with chi-square minimum
# all is weighted with uncertainties

for ioff in range(1, 2001):

    offset = 0.01*ioff

    for j in range(nlines_NIMROD + nlines_SANS2D):

        if q_NIMROD[0] < q_tot_1[j] < max_q_SANS2D:

            chi_square_1[ioff-1] = chi_square_1[ioff-1]+ (f_NIMROD(q_tot_1[j])-offset*f_SANS2D_1(q_tot_1[j]))**2/(f_NIMROD(q_tot_1[j]))

off_min_1 = 10.0
chi_ini_1 = chi_square_1[1999]
offset_list_1 = np.zeros(2000)

for ibis in range(2000):

    offset_list_1[ibis] = 0.01*ibis

    if chi_square_1[ibis] < chi_ini_1:

        chi_ini_1 = chi_square_1[ibis]
        off_min_1 = 0.01*ibis



# Then construction of the merged data by weighted and offset addition of all datasets

q_tot = q_tot_1

dcs_tot = np.zeros(nlines_NIMROD + nlines_SANS2D)
err_dcs_tot = np.zeros(nlines_NIMROD + nlines_SANS2D)

for j in range(nlines_NIMROD + nlines_SANS2D):

    if q_tot_1[j] < q_NIMROD[0]:

        dcs_tot[j] = ((off_min_1 * f_SANS2D_1(q_tot_1[j])/(f_SANS2D_1_err(q_tot_1[j])**2)))/((1/(f_SANS2D_1_err(q_tot_1[j])**2)))
        err_dcs_tot[j] = 1 / np.sqrt(1 / (f_SANS2D_1_err(q_tot_1[j]) ** 2) )

    elif q_tot_1[j] > max_q_SANS2D:

        dcs_tot[j] = f_NIMROD(q_tot_1[j])
        err_dcs_tot[j] = f_NIMROD_err(q_tot_1[j])

    elif max_q_SANS2D >= q_tot_1[j]:

        dcs_tot[j] = (f_NIMROD(q_tot_1[j]) / (f_NIMROD_err(q_tot_1[j]) ** 2) + off_min_1 * f_SANS2D_1(q_tot_1[j]) / (f_SANS2D_1_err(q_tot_1[j]) ** 2)) / (1 / (f_NIMROD_err(q_tot_1[j]) ** 2) + 1 / (f_SANS2D_1_err(q_tot_1[j]) ** 2))
        err_dcs_tot[j] = 1 / np.sqrt(1 / (f_NIMROD_err(q_tot_1[j]) ** 2) + 1 / (f_SANS2D_1_err(q_tot_1[j]) ** 2))

    else:

        dcs_tot[j] = (f_NIMROD(q_tot_1[j]) / (f_NIMROD_err(q_tot_1[j]) ** 2) + off_min_1 * f_SANS2D_1(q_tot_1[j]) / (f_SANS2D_1_err(q_tot_1[j]) ** 2) / (1 / (f_NIMROD_err(q_tot_1[j]) ** 2) + 1 / (f_SANS2D_1_err(q_tot_1[j]) ** 2)+ 1 ))
        err_dcs_tot[j] = 1 / np.sqrt(1 / (f_NIMROD_err(q_tot_1[j]) ** 2) + 1 / (f_SANS2D_1_err(q_tot_1[j]) ** 2)+ 1 )

plt.scatter(q_tot_1, dcs_tot+0.094, c='black', s=2, label='average')
plt.errorbar(q_tot_1, dcs_tot+0.094, err_dcs_tot, c='black', fmt='o', markersize=0.3)

    

# plotting of the Merging data

plt.plot(q_SANS2D, off_min_1 * dcs_SANS2D + 0.094, c='red',label='SANS2D')
plt.plot(q_NIMROD, dcs_NIMROD+0.094, c='blue', label='NIMROD')

plt.xscale('log')
plt.xlim([1e-3, 1])
plt.ylim([1e-2, 1e5])
plt.yscale('log')
plt.xlabel('Q ($\AA^{-1}$)', fontsize=12)
plt.ylabel('Intensity (cm$^{-1}$)', fontsize=12)
plt.legend(frameon=False)
plt.savefig('../Plots/MERGED__witherrors.png', dpi=800)
plt.show()
plt.close()

# export of the Merged data

file_export_merged_alpha = open('../merged_intensity.dat', 'w')

f_merged_alpha = interp1d(q_tot, dcs_tot)
f_merged_alpha_err = interp1d(q_tot, err_dcs_tot)

q_min = q_SANS2D[0]
q_max = 0.7
q_step = 1e-5

line_min = "{:7.5f}".format(q_min) + '\t' + "{:7.5f}".format(f_merged_alpha(q_min)) + '\t' + "{:7.5f}".format(f_merged_alpha_err(q_min) / 100) + ' \n'
file_export_merged_alpha.write(line_min)

for n_q in range(100000):

    q_now = q_min + n_q * q_step
    q_after = q_min + (n_q + 1) * q_step

    for ilines in range(nlines_NIMROD + nlines_SANS2D):

        if q_now < q_tot[ilines] < q_after:
            line_to_write = "{:7.5f}".format(q_tot[ilines]) + '\t' + "{:7.5f}".format(f_merged_alpha(q_tot[ilines])) + '\t' + "{:7.5f}".format(abs(f_merged_alpha_err(q_tot[ilines]))) + '   \n'
            file_export_merged_alpha.write(line_to_write)

file_export_merged_alpha.close()

q_back_min = 0.7
q_back_max = 1.0

npoints = 250
list_back = np.zeros(npoints)

for ipoint in range(npoints):

    list_back[ipoint] = f_merged_alpha(q_back_min+ipoint*0.001)

background = np.mean(list_back)
error_back = np.std(list_back)

print(background, error_back)

plt.scatter(q_tot_1, dcs_tot+0.094, c='black', s=1, label='average')
plt.errorbar(q_tot_1, dcs_tot+0.094, err_dcs_tot, c='black', fmt='o', markersize=0.3)
plt.axhline(y=background+0.094, c='red', label='background = {0}'.format(background))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Q ($\AA^{-1}$)')
plt.ylabel('Intensity (cm$^{-1}$)')
plt.legend(frameon=False)
#plt.savefig('../Plots/Background_Merged_intensity.png', dpi=800)
plt.show()
plt.close()



