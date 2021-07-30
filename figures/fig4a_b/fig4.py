# Setting up of the system to produce fig4.
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

#Set the full path to the flow executable
flowpath='~/kassa-et-al-2021/opm-simulators/build-cmake/bin/flow'

#Import python dependencies
import numpy as np
import meshio
import os
import matplotlib
from matplotlib import pyplot as plt

#Set the parameters for the simulations
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
q=1e-7              #Injection rate [m^3/s]
phi=0.1             #Porosity [-]
tch=1e7             #Characterictic time [s]
b1=1e9              #Dynamic parameter (capillary pressure) [-]
b2=1.8              #Dynamic parameter (capillary pressure) [-]
n1=4.999e5          #Dynamic parameter (relative permeability) [-]
n2=50.              #Dynamic parameter (relative permeability) [-]
nx=640              #Number of cells [-]
dx=2.5              #Length of a grid cell [m]
T=365               #Total simulation time [d]
dt=1.               #Time step to print the results [d]
C=1e-5              #Pore-scale parameter [-]

#Set the three times to print the simulations [d]
tp=[31,181,366]

#Define conversion variables
milli=1e-3          #[-]
darcy=9.8692e-13    #[m^2]
cp=1e-3             #[Pa s]
day=86400           #[s]
atm=101325          #[Pa]

#Create variables used in the simulations
bash = ['1']*2
Sni = []
Snd = []
Wad = []
L=nx*dx
X=np.linspace(dx/2, L-dx/2, nx)
Sw=np.linspace(Srw, 1., 20)
Sew=(Sw-Srw)/(1-Srw)

#Define the saturation functions
def pci(Sw,c):
    return c*Sw**(-1/lambdaa)
def krwi(Sw,E):
     return E*Sw**Lambdaa/(1-Sw+E*Sw**Lambdaa)
def krni(Sw,E):
    return (1-Sw)/(1-Sw+E*Sw**Lambdaa)

#Delete previous simulation files
os.system('rm -r vtk & wait')
os.system('rm -r WACASES & wait')

#Read the input file
a_file = open("WACASE.DATA", "r")
list_of_lines = a_file.readlines()
a_file.close()

#Update the lines using the current input values
list_of_lines[3] = "%d 1 1 /\n" % nx
list_of_lines[51] = "%d*%f /\n" % (nx,dx)
list_of_lines[53] = "%d*1 /\n" % nx
list_of_lines[55] = "%d*1 /\n" % nx
list_of_lines[58] = "%d*160 /\n" % nx
list_of_lines[61] = "%d*%f /\n" % (nx,phi)
list_of_lines[64] = "%d*%f /\n" %(nx,(K/(milli*darcy)))
list_of_lines[67] = "%d*%f /\n" % (nx,(K/(milli*darcy)))
list_of_lines[70] = "%d*%f /\n" % (nx,(K/(milli*darcy)))
for jj in range(20):
    list_of_lines[76+jj] = "%2.4f %2.4f %2.4f %E \n" % (Sw[jj], krwi(Sew[jj],Ei), max(krni(Sew[jj],Ei),0), min(pci(Sew[jj],ci)/atm,100))
list_of_lines[99] = "277.0  1  0.0  %E  0.0 /\n" % (muw/cp)
list_of_lines[102] = "277.0  1  0.0  %E  0.0 /\n" % (mun/cp)
list_of_lines[110] = "%f %f 0 /\n" % (rhon, rhow)
list_of_lines[115] = "%d*1200 /\n" % nx
list_of_lines[119] = "%d*1 /\n" % nx
list_of_lines[123] = "%d*0 /\n" % nx
list_of_lines[133] = "'PROD01' 'PROD' %d 1 1* 'WATER' 0.15/\n" % nx
list_of_lines[138] = "'PROD01' %d 1 1 1 'OPEN' 1* 1*/\n" % nx
list_of_lines[147] = "'INJE01' 'WATER' 'OPEN' 'RATE' %E  1* 1E8 /\n" % (q*day)
list_of_lines[157] = "'INJE01' 'OIL' 'OPEN' 'RATE' %E  1* 1E8 /\n" % (q*day)
list_of_lines[163] = "%d*%f /\n" % (int(T/dt),dt)

#Set the simulation
os.system('mkdir vtk & wait')
os.system('mkdir WACASES & wait')
a_file = open("WACASES/WACASE-0.DATA", "w")
a_file.writelines(list_of_lines)
a_file.close()
os.system('mkdir vtk/vtk-0 & wait')
bash[0]="%s WACASES/WACASE-0.DATA --output-dir=vtk/vtk-0 --enable-vtk-output=true --enable-ecl-output=false  --initial-time-step-in-days=0.0001 --solver-max-restarts=20 --solver-max-time-step-in-days=.1 --enable-wa=true --beta=%E --eta=%E --ei=%f --ef=%f --ci=%E --cf=%E --lambda=%f --llambda=%f --srw=%f" % (flowpath,1e16,0,Ei,Ef,ci,ci,lambdaa,Lambdaa,Srw)
os.system('mkdir vtk/vtk-1 & wait')
bash[1]="%s WACASES/WACASE-0.DATA --output-dir=vtk/vtk-1 --enable-vtk-output=true --enable-ecl-output=false  --initial-time-step-in-days=0.0001 --solver-max-restarts=20 --solver-max-time-step-in-days=.1 --enable-wa=true --beta=%E --eta=%E --ei=%f --ef=%f --ci=%E --cf=%E --lambda=%f --llambda=%f --srw=%f" % (flowpath,b1*C**b2,-n1*C+n2,Ei,Ef,ci,cf,lambdaa,Lambdaa,Srw)

#Create the .bash file and run the simulation
a_file = open("fig4.bash", "w")
a_file.writelines(bash[0]+" & "+ bash[1]+ " & wait\n")
a_file.close()
os.system("chmod u+x ./fig4.bash")
os.system('./fig4.bash')

#Obtain the final Sn and Wa from the vtk output files
for i in tp:
    mesh = meshio.read("vtk/vtk-0/WACASE-0-%05d.vtu" % i)
    for row in mesh.cell_data['saturation_oil']:
        sni=row
    Sni.append(sni)
    mesh = meshio.read("vtk/vtk-1/WACASE-0-%05d.vtu" % i)
    for row in mesh.cell_data['saturation_oil']:
        snd=row
    Snd.append(snd)
    for row in mesh.cell_data['wa']:
        wad=row
    Wad.append(wad)

#Plot the results and save them as fig4a.eps
lw=1
plt.figure(figsize=(6., 5.5), dpi=80)
plt.rc('font', size=12)
axes=plt.subplot(1, 1, 1)
plt.plot(X, Wad[0][:], color=[0,0,0], linewidth=lw, linestyle="-",label=r"Dynamic k$_r$ P$_c$, t=30 d")
plt.plot(X, Wad[1][:], color=[0,0,0], linewidth=lw, linestyle=":",label=r"Dynamic k$_r$ P$_c$, t=180 d")
plt.plot(X, Wad[2][:], color=[0,0,0], linewidth=lw, linestyle="--",label=r"Dynamic k$_r$ P$_c$, t=365 d")
plt.yscale('log')
plt.xlim([0,L])
plt.ylim([1e-8,1])
plt.xlabel('x [m]')
plt.ylabel('$\overline{\chi}$ [-]')
plt.title('(a)')
plt.grid()
plt.legend(loc='best', prop={'size': 9})
plt.savefig('fig4a.eps', format='eps')
plt.show()

#Plot the results and save them as fig4b.eps
lw=1
plt.figure(figsize=(6., 5.5), dpi=80)
plt.rc('font', size=12)
axes=plt.subplot(1, 1, 1)
plt.plot(X, Sni[0][:], color=[1,.6,.6], linewidth=lw, linestyle="-",label="Static k$_r$ P$_c$, t=30 d")
plt.plot(X, Sni[1][:], color=[1,.6,.6], linewidth=lw, linestyle=":",label="Static k$_r$ P$_c$, t=180 d")
plt.plot(X, Sni[2][:], color=[1,.6,.6], linewidth=lw, linestyle="--",label="Static k$_r$ P$_c$, t=365 d")
plt.plot(X, Snd[0][:], color=[0,0,0], linewidth=lw, linestyle="-",label=r"Dynamic k$_r$ P$_c$, t=30 d")
plt.plot(X, Snd[1][:], color=[0,0,0], linewidth=lw, linestyle=":",label=r"Dynamic k$_r$ P$_c$, t=180 d")
plt.plot(X, Snd[2][:], color=[0,0,0], linewidth=lw, linestyle="--",label=r"Dynamic k$_r$ P$_c$, t=365 d")
plt.yscale('log')
plt.xlim([0,L])
plt.ylim([1e-8,1])
plt.xlabel('x [m]')
plt.ylabel('S$_n$ [-]')
plt.title('(b)')
plt.grid()
plt.legend(loc='upper right', prop={'size': 9})
plt.savefig('fig4b.eps', format='eps')
plt.show()
