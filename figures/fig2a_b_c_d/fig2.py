# Setting up of the system to produce fig2a, fig2b, fig2c, and fig2d.
#{
#Copyright 2021, NORCE Norwegian Research Centre AS, Computational
#Geosciences and Modeling.

#This file is part of the ad-wa module.

#ad-wa is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#ad-wa is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this file.  If not, see <http://www.gnu.org/licenses/>.
#}

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

#Set the parameters for the saturation functions
ci      = 1e4       #Initial entry pressure [Pa]
cf      = 1e2       #Final entry pressure [Pa]
lambdaa = 3.6       #Fitting parameter (relative permeability) [-]
Lambdaa = 1.3       #Fitting parameter (capillary pressure) [-]
Ei      = 0.48      #Initial wetting-parameter (relative permeability) [-]
Ef      = 3.37      #Final wetting-parameter (relative permeability) [-]
rhon    = 716.7     #Non-wetting density [kg/m^3]
mun     = 5.916e-5  #Non-wetting viscocity [Pa s]
rhow    = 1050      #Wetting density [kg/m^3]
muw     = 6.922e-4  #Wetting viscoscity [Pa s]
K       = 1e-10     #Intrinsic permeability [m^2]
q       = 1e-3      #Injection rate [m^3/s]
phi     = 0.1       #Porosity [-]
g       = 9.81      #Gravity [m/s^2]

#Define the saturation functions
def pc(Sew, c):
    return c * Sew ** (-1 / lambdaa)
def krw(Sew, E):
    return E * Sew ** Lambdaa / ((1 - Sew) + E * Sew ** Lambdaa)
def krn(Sew, E):
    return (1 - Sew) / ((1 - Sew) + E * Sew ** Lambdaa)
def fn(Sew, E, v):
    return (1 - (K / q) * (krw(Sew, E) / muw) * (rhon - rhow) * g * np.sin(v)) \
           / (1 + ((krw(Sew, E) / muw) / (krn(Sew, E) / mun)))

#Create variables used in the plotting
Sew   = np.arange(0.0001, 1 - 0.0001, 0.0001)
pcv   = np.vectorize(pc)
krwv  = np.vectorize(krw)
krnv  = np.vectorize(krn)
fnv   = np.vectorize(fn)

#Plot the curves and save them as fig2a.eps
lw = 4
plt.figure(figsize = (6, 5.5), dpi = 80)
plt.rc('font', size = 12)
axes = plt.subplot(1, 1, 1)
plt.plot(Sew, pcv(Sew, ci / 1e6), color = [0, 0, 0], linewidth = lw,
         linestyle = "--", label = "Initial-wetting function")
plt.plot(Sew, pcv(Sew, cf / 1e6), color = [1, 0, 0],
         linewidth = lw, linestyle = "-", label = "Final-wetting function")
plt.xlim([0, 1])
plt.ylim([1e-4, 1e-1])
plt.yscale('log')
plt.xlabel('S$_{ew}$ [-]')
plt.ylabel('P$_c$ [MPa]')
plt.title('(a)')
plt.grid()
plt.legend(loc = 'best')
plt.savefig('fig2a.eps', format = 'eps')

#Plot the curves and save them as fig2b.eps
plt.figure(figsize = (6, 5.5), dpi = 80)
plt.rc('font', size = 12)
axes = plt.subplot(1, 1, 1)
plt.plot(Sew, krwv(Sew, Ei), color = [0, 0, 0], linewidth = lw,
         linestyle = "-", label = "Initial-wetting function")
plt.plot(Sew, krnv(Sew, Ei), color = [0, 0, 0], linewidth = lw,
         linestyle = "--",label = "Initial-wetting function")
plt.plot(Sew, krwv(Sew, Ef), color = [1, 0, 0], linewidth = lw,
         linestyle = "-", label = "Final-wetting function")
plt.plot(Sew, krnv(Sew, Ef), color = [1, 0, 0], linewidth = lw,
         linestyle = ":", label = "Final-wetting function")
plt.annotate('', xy = (0.76, 0.6), xytext=(0.7, 0.7), arrowprops = dict(
             facecolor = 'blue', shrink = 0.01, width = 0.5, headwidth = 7.5))
plt.annotate('', xy = (0.6, 0.8), xytext = (0.7, 0.7), arrowprops = dict(
             facecolor = 'blue', shrink = 0.01, width = 0.5, headwidth = 7.5))
plt.annotate('', xy = (0.76, 0.4), xytext = (0.7, 0.3), arrowprops = dict(
             facecolor = 'yellow', shrink = 0.01, width = 0.5, headwidth = 7.5))
plt.annotate('', xy = (0.6, 0.2), xytext = (0.7, 0.3), arrowprops = dict(
             facecolor = 'yellow', shrink = 0.01, width = 0.5, headwidth = 7.5))
plt.text(0.7, 0.03, 'Non-wetting phase', rotation = -42, fontsize = 12)
plt.text(0.7, 0.7, 'Wetting phase', rotation = 50, fontsize = 12)
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('S$_{ew}$ [-]')
plt.ylabel(r'k$_{r \alpha}$ [-]')
plt.title('(b)')
plt.grid()
plt.legend(loc = 'best')
plt.savefig('fig2b.eps', format = 'eps')

#Plot the curves and save them as fig2c.eps
plt.figure(figsize=(6, 5.5), dpi = 80)
plt.rc('font', size = 12)
axes = plt.subplot(1, 1, 1)
plt.plot(Sew, fn(Sew, Ei, 0), color = [0, 0, 0], linewidth = lw,
         linestyle = "--", label = "Initial-wetting function")
plt.plot(Sew, fn(Sew, Ef, 0), color = [1, 0, 0], linewidth = lw,
         linestyle = "-", label = "Final-wetting function")
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('S$_{ew}$ [-]')
plt.ylabel('f$_n$ [-]')
plt.title('(c)')
plt.grid()
plt.legend(loc = 'best')
plt.savefig('fig2c.eps', format = 'eps')

#Plot the curves and save them as fig2d.eps
plt.figure(figsize = (6, 5.5), dpi = 80)
plt.rc('font', size = 12)
axes = plt.subplot(1, 1, 1)
plt.plot(Sew, fn(Sew, Ei, np.pi / 2), color = [0, 0, 0], linewidth = lw,
         linestyle = "--", label = "Initial-wetting function")
plt.plot(Sew, fn(Sew, Ef, np.pi / 2), color = [1, 0, 0], linewidth = lw,
         linestyle = "-", label = "Final-wetting function")
plt.xlim([0, 1])
plt.ylim([0, 1.2])
plt.xlabel('S$_{ew}$ [-]')
plt.ylabel('f$_n$ [-]')
plt.title('(d)')
plt.grid()
plt.legend(loc = 'best')
plt.savefig('fig2d.eps', format = 'eps')
