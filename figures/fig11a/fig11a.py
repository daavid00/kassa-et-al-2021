# Setting up of the system to produce fig11a.
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
from scipy.optimize import curve_fit

#Set the parameters for the scaling
cic=1e4             #Reference entry pressure (caprock) [Pa]
Kc=1e-16            #Reference permeability (caprock) [m^2]
C=1e-5              #Reference pore-scale parameter [-]
tref=10000.         #Reference time [y]
Mmax=2733.          #Maximum CO2 mass [kg/m^2]

#Define the cases for the scaling
M=3
ii=0
b=[1.,5.,10.]
c=[1.,0.75,0.5]
d=[1.,2.5,5.]
H = []
for j in range(M):
    for i in range(M):
        a = []
        a.append(b[0]*cic)
        a.append(c[i]*Kc)
        a.append(d[j]*C)
        H.append(a)
        ii=ii+1
    for i in range(M):
        a = []
        a.append(b[i]*cic)
        a.append(c[0]*Kc)
        a.append(d[j]*C)
        H.append(a)
        ii=ii+1
P = []
for item in H:
    if item not in P:
             P.append(item)
N = len(P)

#Create variables used in the scaling
CO2 = []
tt = []
Tp = []
th = [0]*N
CO2h = [0]*N
tth = [0]*N
CO2scaled = [0]*N
CO2ref = []
ttref = []
xci = []
xco2 = []
xkc = []
xc = []
yco2 = []
tco2 = []
mco2 = []

#First run the simulations for fig10a, then copy the produced folder 'data-for-fig11a' here
CO2ref.append(np.loadtxt('data-for-fig11a/wa-fig10a-co2-00000.csv', delimiter=','))
ttref.append(np.loadtxt('data-for-fig11a/wa-fig10a-time-00000.csv', delimiter=','))
for i in range(N):
    CO2.append(np.loadtxt('data-for-fig11a/wa-fig10a-co2-%05d.csv' % (i+1), delimiter=','))
    tt.append(np.loadtxt('data-for-fig11a/wa-fig10a-time-%05d.csv' % (i+1), delimiter=','))

#Define the scaling functions
def timehat(X, a0, a1):
    xci,xc = X
    return a0*(xci/cic)*(xc/C)**a1
def co2hat(X, a2, a3, a4, a5, a6):
    xci,xkc,xc,xco2 = X
    return ((xkc/Kc)**a2*(xc/C)**a3+a4*((xci)/cic)**a5*(xc/C)**a6)*xco2
def modelhat(X, psi0):
    t = X
    return psi0*t

#Estimate the parameters for the time scaling
txci=[row[0] for row in P]
txc=[row[2] for row in P]
for i in range(N):
    th[i]=(len(CO2[i][tuple([CO2[i][:]<1e-10])]))/(tref)
p0=1e5 ,0.1
parameters, covariance = curve_fit(timehat, (txci,txc), th, p0)
a0 = parameters[0]
a1 = parameters[1]

#Scale the time
for i in range(N):
    CO2h[i]=CO2[i][tuple([tt[i][:]/tref-timehat((P[i][0],P[i][2]),a0,a1)>0])]
    tth[i]=np.linspace(0, len(CO2h[i][:])/(tref), len(CO2h[i][:]))

#Estimate the parameters for the CO2 scaling
for i in range(N):
    for j in range(len(CO2h[i][:])):
        xci.append(P[i][0])
        xkc.append(P[i][1])
        xc.append(P[i][2])
        xco2.append(CO2h[i][j]/Mmax)
        yco2.append(CO2ref[0][j]/Mmax)
p0=-0.7, 0.08, 0.01, 1.0, 0.5
parameters, covariance = curve_fit(co2hat, (xci,xkc,xc,xco2), yco2, p0)
a2 = parameters[0]
a3 = parameters[1]
a4 = parameters[2]
a5 = parameters[3]
a6 = parameters[4]

#Scale the CO2
for i in range(N):
    CO2scaled[i]=co2hat((P[i][0],P[i][1],P[i][2],1.0),a2,a3,a4,a5,a6)*(np.array(CO2h[i][:],dtype=float)/Mmax)

#Estimate the parameters for the model
parameters, covariance = curve_fit(modelhat, np.array(ttref[0][:],dtype=float)/tref, (CO2ref[0][:])/Mmax)
psi0 = parameters[0]

#Plot and save the results for fig11a.eps
lw=5
nn=5
plt.figure(figsize=(8, 6), dpi=80)
plt.rc('font', size=12)
plt.rc('legend', fontsize=9)
axes=plt.subplot(1, 1, 1)

plt.plot(tth[0], CO2scaled[0], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="o",markevery=nn, label=r"K=$10^{-16}$ m$^2$, c$^i$=$10^4$ Pa, C=$10^{-5}$")
plt.plot(tth[1], CO2scaled[1], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="s",markevery=nn, label=r"K=$7.5\times 10^{-17}$ m$^2$")
plt.plot(tth[2], CO2scaled[2], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="d",markevery=nn, label=r"K=$5\times 10^{-17}$ m$^2$")
plt.plot(tth[3], CO2scaled[3], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="P",markevery=nn, label=r"$c^i=5\times 10^4$ Pa")
plt.plot(tth[4], CO2scaled[4], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="v",markevery=nn, label=r"$c^i=10^5$ Pa")
plt.plot(tth[5], CO2scaled[5], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="o",markevery=nn, label=r"K=$10^{-16}$ m$^2$, c$^i$=$10^4$ Pa, C=$2.5\times 10^{-5}$")
plt.plot(tth[6], CO2scaled[6], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="s",markevery=nn, label=r"K=$7.5\times 10^{-17}$ m$^2$")
plt.plot(tth[7], CO2scaled[7], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="d",markevery=nn, label=r"K=$5\times 10^{-17}$ m$^2$")
plt.plot(tth[8], CO2scaled[8], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="P",markevery=nn, label=r"$c^i=5\times 10^4$ Pa")
plt.plot(tth[9], CO2scaled[9], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="v",markevery=nn, label=r"$c^i=10^5$ Pa")
plt.plot(tth[10], CO2scaled[10], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="o",markevery=nn, label=r"K=$10^{-16}$ m$^2$, c$^i$=$10^4$ Pa, C=$5\times 10^{-5}$")
plt.plot(tth[11], CO2scaled[11], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="s",markevery=nn, label=r"K=$7.5\times 10^{-17}$ m$^2$")
plt.plot(tth[12], CO2scaled[12], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="d",markevery=nn, label=r"K=$5\times 10^{-17}$ m$^2$")
plt.plot(tth[13], CO2scaled[13], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="P",markevery=nn, label=r"$c^i=5\times 10^4$ Pa")
plt.plot(tth[14], CO2scaled[14], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="v",markevery=nn, label=r"$c^i=10^5$ Pa")
plt.plot(np.array(ttref[0][:],dtype=float)/tref, modelhat(np.array(ttref[0][:],dtype=float)/tref,psi0), linewidth=lw, linestyle="--", color=[0.5,0.5,0.5], label=r"Model")

plt.xlim([0,0.01])
plt.ylim([0,.15])
plt.xlabel('$\hat{t}$ [-]')
plt.ylabel('$\hat{M}_{CO_{2}}$ [-]')
plt.grid()
matplotlib.pyplot.grid(True, which="both")
plt.legend(loc='upper left')
plt.savefig('fig11a.eps', format='eps')
plt.show()
