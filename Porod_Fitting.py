# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from matplotlib import cm
import pandas as pd
import glob
import matplotlib as mpl
import matplotlib.pyplot as mplplot
from sklearn import preprocessing
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit


###################### SETTING UP DATA

data = pd.read_csv("C:/Users/jch13924/Documents/PhD/Neutron/NIMROD/cycle_15_1/Sample_4_30K/isothermal_steps/mints/Dep30K_heating_30K.mint01",
                   sep='\s+',
                   skiprows=14,
                   nrows=126,
                   names=["Q", "DCS", "DCS err"],
                   )

#print(data)

I_data = data["DCS"] * 0.094 * 1e-4
Ierr_data = data["DCS err"] * 0.094 * 1e-4

q_data = data["Q"]

y_data = I_data*q_data**4
yerr_data = Ierr_data*q_data**4


porod_data = pd.read_csv("C:/Users/jch13924/Documents/PhD/Neutron/NIMROD/cycle_15_1/Sample_4_30K/isothermal_steps/mints/Dep30K_heating_30K.mint01",
                   sep='\s+',
                   skiprows=107,
                   nrows= 18,
                   names=["Q", "DCS", "DCS err"],
                   )


#print(porod_data)

I = porod_data["DCS"] * 0.094 * 1e-4
Ierr = porod_data["DCS err"] * 0.094 * 1e-4

q = porod_data["Q"]

y = I*q**4
yerr = Ierr*q**4
 

###################### FITTING

weight = yerr / y
p0 = 1e-8

result, covariant = curve_fit(lambda x, p0: p0, q, y, p0=p0, sigma=weight, absolute_sigma=True, maxfev=10000)

result_error = np.sqrt(((y - result)**2).sum() / (len(y)-1))

SSA = result / (2.0*np.pi*6.38e-6**2.0)
SSA_err = result_error / (2.0*np.pi*6.38e-6**2.0)

K_line = np.array(np.zeros(len(y)) + result)


###################### PLOTTING


plt.figure(dpi=1200)
plt.figure(figsize= (8,6))
mpl.rcParams["axes.linewidth"]	= 1.5
fig,ax = plt.subplots()

ax.plot(q_data, y_data, color='darkred')
ax.fill_between(q_data, y_data-yerr_data, y_data+yerr_data, color='darkred', alpha=0.5, zorder=2)

ax.plot(q, y, color='k', linestyle='--')
ax.plot(q, K_line, 'k-')
ax.plot(q, K_line+result_error,'k--')
ax.plot(q, K_line-result_error, 'k--')

ax.text(0.7, 0.08, "K = %.1e $\pm$ %.1e" %(result, result_error), transform=ax.transAxes, horizontalalignment='left', verticalalignment='bottom', color="k", fontsize=9)
ax.text(0.7, 0.04, "SSA = %.2f $\pm$ %.2f" %(SSA, SSA_err), transform=ax.transAxes, horizontalalignment='left', verticalalignment='bottom', color="k", fontsize=9)


ax.tick_params(which="major", labelsize=12, width=1.5, length=6)
ax.tick_params(which="minor", labelsize=12, width=1.5, length=4)
ax.set_xlabel('Q ($\AA$$^-$$^1$)', fontsize=12)
ax.set_ylabel('I(Q)*Q4', fontsize=12)

plt.show()









