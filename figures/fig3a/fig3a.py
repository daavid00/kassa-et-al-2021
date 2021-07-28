# Setting up of the system to produce fig3a.
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
nx=500              #Number of cells [-]
dx=2.5              #Length of a grid cell [m]
T=365               #Total simulation time [d]
dt=73.              #Time step to print the results [d]

#Define conversion variables
milli=1e-3          #[-]
darcy=9.8692e-13    #[m^2]
cp=1e-3             #[Pa s]
day=86400           #[s]
atm=101325          #[Pa]

#Define the cases for the static states
N=4
P = []
a = []
a.append(Ei)
a.append(ci)
P.append(a)
a = []
a.append(Ef)
a.append(ci)
P.append(a)
a = []
a.append(Ei)
a.append(cf)
P.append(a)
a = []
a.append(Ef)
a.append(cf)
P.append(a)

#Create variables used in the simulations
Sn = []
bashi = ['1']*N
Oki = [0]*N
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
list_of_lines[99] = "277.0  1  0.0  %E  0.0 /\n" % (muw/cp)
list_of_lines[102] = "277.0  1  0.0  %E  0.0 /\n" % (mun/cp)
list_of_lines[110] = "%f %f 0 /\n" % (rhon, rhow)
list_of_lines[115] = "%d*%d /\n" % (nx,nx*2.)
list_of_lines[119] = "%d*1 /\n" % nx
list_of_lines[123] = "%d*0 /\n" % nx
list_of_lines[133] = "'PROD01' 'PROD' %d 1 1* 'WATER' 0.15/\n" % nx
list_of_lines[138] = "'PROD01' %d 1 1 1 'OPEN' 1* 1*/\n" % nx
list_of_lines[147] = "'INJE01' 'WATER' 'OPEN' 'RATE' %E  1* 1E8 /\n" % (q*day)
list_of_lines[157] = "'INJE01' 'OIL' 'OPEN' 'RATE' %E  1* 1E8 /\n" % (q*day)
list_of_lines[163] = "%d*%f /\n" % (int(T/dt),dt)

#Set the simulations for the different static cases
os.system('mkdir vtk & wait')
os.system('mkdir WACASES & wait')
for i in range(N):
    for jj in range(20):
        list_of_lines[76+jj] = "%2.4f %2.4f %2.4f %E \n" % (Sw[jj], krwi(Sew[jj],P[i][0]), max(krni(Sew[jj],P[i][0]),0), min(pci(Sew[jj],P[i][1])/atm,100))
    a_file = open("WACASES/WACASE-%05d.DATA" % i, "w")
    a_file.writelines(list_of_lines)
    a_file.close()
    os.system('mkdir vtk/vtk-%05d & wait' % i)
    #Set the values for the simulation
    bashi[i]="%s WACASES/WACASE-%05d.DATA --output-dir=vtk/vtk-%05d --enable-vtk-output=true --enable-ecl-output=false --initial-time-step-in-days=0.0001 --solver-max-restarts=20 --solver-max-time-step-in-days=1. --enable-wa=false" % (flowpath,i,i)

#Create the .bash file and run the simulations
bash=bashi[0]+" & "+bashi[1]+" & "+bashi[2]+" & "+bashi[3]+" & wait\n"
a_file = open("fig3a.bash", "w")
a_file.writelines(bash)
a_file.close()
os.system("chmod u+x ./fig3a.bash")
os.system('./fig3a.bash')

#Obtain the final Sn from the vtk output files
for i in range(N):
    a_file = open("vtk/vtk-%05d/WACASE-%05d.pvd" % (i,i),"r")
    list_of_lines = a_file.readlines()
    ss=list_of_lines[-3]
    Oki[i]=int(ss[-13:-8])
    mesh = meshio.read("vtk/vtk-%05d/WACASE-%05d-%05d.vtu" % (i,i,int(ss[-13:-8])))
    for row in mesh.cell_data['saturation_oil']:
        snd=row
    Sn.append(snd)

#Check that all different simulations finished (the number of output files should be the same)
print(Oki)

#Plot the results and save them as fig3a.eps
lw=4
plt.figure(figsize=(5.5, 5.5), dpi=80)
plt.rc('font', size=12)
axes=plt.subplot(1, 1, 1)
plt.plot(X, Sn[0][:], color=[1,0,.5], linewidth=lw, linestyle="-", label="Initial-wetting k$_r$, Initial-wetting P$_c$")
plt.plot(X, Sn[1][:], color=[0,.5,1], linewidth=lw, linestyle="-", label="Final-wetting k$_r$, Initial-wetting P$_c$")
plt.plot(X, Sn[2][:], color=[.8,.8,0], linewidth=lw, linestyle="-", label="Initial-wetting k$_r$, Final-wetting P$_c$")
plt.plot(X, Sn[3][:], color=[.3,0,.6], linewidth=lw, linestyle="-", label="Final-wetting k$_r$, Final-wetting P$_c$")
plt.xlim([0,L])
plt.ylim([0,0.425])
plt.xlabel('x [m]')
plt.ylabel('S$_n$ [-]')
plt.grid()
plt.title('(a)')
plt.legend(loc='best')
plt.savefig('fig3a.eps', format='eps')
plt.show()
