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
ci=1e4              #Initial entry pressure [Pa]
cf=1e2              #Final entry pressure [Pa]
lambdaa=3.6         #Fitting parameter (relative permeability) [-]
Lambdaa=1.3         #Fitting parameter (capillary pressure) [-]
Ei=0.48             #Initial wetting-parameter (relative permeability) [-]
Ef=3.37             #Final wetting-parameter (relative permeability) [-]
rhon=716.7          #Non-wetting density [kg/m^3]
mun=5.916e-5        #Non-wetting viscocity [Pa s]
rhow=1050           #Wetting density [kg/m^3]
muw=6.922e-4        #Wetting viscoscity [Pa s]
K=1e-10             #Intrinsic permeability [m^2]
Srw=0.2             #Residual wetting saturation [-]
Srn=0.0             #Residual non-wetting saturation [-]
q=1e-3              #Injection rate [m^3/s]
phi=0.1             #Porosity [-]
g=9.81              #Gravity [m/s^2]

#Define the saturation functions
def pc(Sn,c):
    return c*(1-Sn)**(-1/lambdaa)
def krw(Sn,E):
    return E*(1-Sn)**Lambdaa/(Sn+E*(1-Sn)**Lambdaa)
def krn(Sn,E):
    return Sn/(Sn+E*(1-Sn)**Lambdaa)
def fn(Sn,E,v):
    return (1-(K/q)*(krw(Sn,E)/muw)*(rhon-rhow)*g*np.sin(v))/(1+((krw(Sn,E)/muw)/(krn(Sn,E)/mun)))

#Create variables used in the plotting
Snmax=1.0-Srw
Sn=np.arange(Srn+0.0001, Snmax-0.0001, 0.0001)
Sne=(Sn-Srn)/(1-Srw-Srn)    
pcv=np.vectorize(pc)
krwv=np.vectorize(krw)
krnv=np.vectorize(krn)
fnv=np.vectorize(fn)

#Plot the curves and save them as fig2a.eps
lw=4
plt.figure(figsize=(6, 5.5), dpi=80)
plt.rc('font', size=12)
axes=plt.subplot(1, 1, 1)
plt.plot(Sn,pcv(Sne,ci/1e6), color=[0,0,0], linewidth=lw, linestyle="--", label="Initial-wetting function")
plt.plot(Sn,pcv(Sne,cf/1e6), color=[1,0,0], linewidth=lw, linestyle="-", label="Final-wetting function")
plt.xlim([Srn,Snmax])
plt.ylim([1e-4,1e-1])
plt.yscale('log')
plt.xlabel('S$_n$ [-]')
plt.ylabel('P$_c$ [MPa]')
plt.title('(a)')
plt.grid()
plt.legend(loc='upper left')
plt.savefig('fig2a.eps', format='eps')

#Plot the curves and save them as fig2b.eps
plt.figure(figsize=(6, 5.5), dpi=80)
plt.rc('font', size=12)
axes=plt.subplot(1, 1, 1)
plt.plot(Sn,krwv(Sne,Ei), color=[0,0,0], linewidth=lw, linestyle="-", label="Initial-wetting function")
plt.plot(Sn,krnv(Sne,Ei), color=[0,0,0], linewidth=lw, linestyle="--", label="Initial-wetting function")
plt.plot(Sn,krwv(Sne,Ef), color=[1,0,0], linewidth=lw, linestyle="-", label="Final-wetting function")
plt.plot(Sn,krnv(Sne,Ef), color=[1,0,0], linewidth=lw, linestyle=":", label="Final-wetting function")
plt.xlim([Srn,Snmax])
plt.ylim([0,1])
plt.xlabel('S$_n$ [-]')
plt.ylabel(r'k$_{r \alpha}$ [-]')
plt.title('(b)')
plt.grid()
plt.legend(loc='best')
plt.savefig('fig2b.eps', format='eps')

#Plot the curves and save them as fig2c.eps
plt.figure(figsize=(6, 5.5), dpi=80)
plt.rc('font', size=12)
axes=plt.subplot(1, 1, 1)
plt.plot(Sn,fn(Sne,Ei,0), color=[0,0,0], linewidth=lw, linestyle="--", label="Initial-wetting function")
plt.plot(Sn,fn(Sne,Ef,0), color=[1,0,0], linewidth=lw, linestyle="-", label="Final-wetting function")
plt.xlim([Srn,Snmax])
plt.ylim([0,1])
plt.xlabel('S$_n$ [-]')
plt.ylabel('f$_n$ [-]')
plt.title('(c)')
plt.grid()
plt.legend(loc='best')
plt.savefig('fig2c.eps', format='eps')

#Plot the curves and save them as fig2d.eps
plt.figure(figsize=(6, 5.5), dpi=80)
plt.rc('font', size=12)
axes=plt.subplot(1, 1, 1)
plt.plot(Sn,fn(Sne,Ei,np.pi/2), color=[0,0,0], linewidth=lw, linestyle="--", label="Initial-wetting function")
plt.plot(Sn,fn(Sne,Ef,np.pi/2), color=[1,0,0], linewidth=lw, linestyle="-", label="Final-wetting function")
plt.xlim([Srn,Snmax])
plt.ylim([0,1.2])
plt.xlabel('S$_n$ [-]')
plt.ylabel('f$_n$ [-]')
plt.title('(d)')
plt.grid()
plt.legend(loc='best')
plt.savefig('fig2d.eps', format='eps')
