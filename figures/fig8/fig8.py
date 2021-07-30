# Setting up of the system to produce fig8.
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

#Set the full path to the wa executable
flowwa='~/kassa-et-al-2021/opm-models/build-cmake/bin/wa'

#Import python dependencies
import numpy as np
import meshio
import os
import matplotlib
from matplotlib import pyplot as plt

#Set the parameters for the simulations
lambdaa=3.6         #Fitting parameter (capillary pressure) [-]
Lambdaa=1.3         #Fitting parameter (relative permeability) [-]
Ei=0.48             #Initial wetting-parameter (relative permeability) [-]
Ef=3.37             #Final wetting-parameter (relative permeability) [-]
Srw=0.2             #Residual wetting saturation [-]
Srn=0.2             #Residual non-wetting saturation [-]
cia=1e2             #Initial entry pressure (aquifer) [Pa]
cic=1e4             #Initial entry pressure (caprock) [Pa]
cfa=1e2             #Final entry pressure (aquifer) [Pa]
cfc=1e2             #Final entry pressure (caprock) [Pa]
rhon=716.7          #Non-wetting density [kg/m^3]
Ka=1e-10            #Intrinsic permeability (aquifer) [m^2]
Kc=1e-16            #Intrinsic permeability (caprock) [m^2]
phia=0.2            #Porosity (aquifer) [-]
phic=0.2            #Porosity (caprock) [-]
nz=100              #Number of cells (z) [-]
dx=10               #Size of a grid cell (x) [m]
dz=1                #Size of a grid cell (z) [m]
h=70.               #Height of the aquifer [m]
S0w=0.5             #Initial wetting saturation (aquifer) [-]
X0n=5e-3            #Initial co2 mole fraction (brine) [-]
T=4000              #Total simulation time [y]
dt=3650.            #Time step to print the results [d]

#Define conversion variables
day=86400           #[s]
year=365*day        #[s]

#Create the file for the printed times
X = ['1']*(int(T*year/(dt*day)))
for i in range(int(T*year/(dt*day))):
    X[i]="%d\n" % ((i+1)*dt*day)
a_file = open("writetimes.DATA", "w")
a_file.writelines(X)
a_file.close()

#Update the file for the grid
a_file = open("wa.dgf", "r")
list_of_lines = a_file.readlines()
a_file.close()
a_file = open("wa.dgf", "w")
list_of_lines[4] = "%d %d %d %% \n" % (dx,1,nz*dz)
list_of_lines[5] = "%d %d %d %% \n" % (1,1,nz)
a_file.writelines(list_of_lines)
a_file.close()

#Define the cases for the simulations
N = 6
ii=0
P = []
a = []
a.append(Ei)
a.append(Kc)
P.append(a)
a = []
a.append(Ef)
a.append(Kc)
P.append(a)
a = []
a.append(Ei)
a.append(0.75*Kc)
P.append(a)
a = []
a.append(Ef)
a.append(0.75*Kc)
P.append(a)
a = []
a.append(Ei)
a.append(0.5*Kc)
P.append(a)
a = []
a.append(Ef)
a.append(0.5*Kc)
P.append(a)

#Create variables used in the simulations
bashi = [0]*N
Oki = [0]*N
CO2 = []

#Delete previous simulation files
os.system('rm -r vtk & wait')

#Set the different simulations
os.system('mkdir vtk')
for i in range(N):
    os.system('mkdir vtk/vtk-%05d & wait' % i)
    bashi[i]="%s --output-dir=vtk/vtk-%05d --wa-vtk-time-steps-file=writetimes.DATA --initial-time-step-size=.01 --max-time-step-size=2764800 --max-time-step-divisions=20  --end-time=%d --enable-wa=true --beta=%E --eta=%E --ci-c=%f --cf-c=%f --ci=%f --cf=%f --ei=%f --ef=%f --lambda=%f --llambda=%f --srw=%f --srn=%f --s0w=%f --x0n=%f --k-c=%E --k=%E --phi-c=%f --phi=%f --fine-layer-bottom=%f " % (flowwa,i,T*year,1e15,0,cfc,cfc,cia,cfa,P[i][0],P[i][0],lambdaa,Lambdaa,Srw,Srn,S0w,X0n,P[i][1],Ka,phic,phia,h)

#Create the .bash file and run the simulations
bash=bashi[0]+" & "+bashi[1]+" & "+bashi[2]+" & "+bashi[3]+" & "+bashi[4]+" & "+bashi[5]+" & wait\n"
a_file = open("fig8.bash", "w")
a_file.writelines(bash)
a_file.close()
os.system("chmod u+x ./fig8.bash")
os.system('./fig8.bash')

#Compute the CO2 in the caprock from the vtk output files
mesh = meshio.read("vtk/vtk-00000/wa_ncp_ecfv-00000.vtu")
for row in mesh.cell_data['saturation_gas']:
    sn=row
inx=np.where(sn<1e-10)
for i in range(N):
    a = []
    a_file = open("vtk/vtk-%05d/wa_ncp_ecfv.pvd" % i,"r")
    list_of_lines = a_file.readlines()
    ss=list_of_lines[-3]
    Oki[i]=int(ss[-13:-8])
    for j in range(int(ss[-13:-8])+1):
        mesh = meshio.read("vtk/vtk-%05d/wa_ncp_ecfv-%05d.vtu" % (i,j))
        for row in mesh.cell_data['saturation_gas']:
            snd=row
        a_file.close()
        a.append(sum(rhon*phic*snd[inx]))
    CO2.append(a)

#Check that all different simulations finished (the number of output files should be the same (+-1))
print(Oki)

#Plot and save the results for fig8.eps
lw=2
plt.figure(figsize=(8, 6), dpi=512)
plt.rc('font', size=12)
plt.rc('legend', fontsize=9)
axes=plt.subplot(1, 1, 1)
plt.plot(np.linspace(0, len(CO2[0][:])*dt*day/year, len(CO2[0][:])), CO2[0][:], color=[.3,0,.6], linewidth=lw, linestyle="-", label=r"Initial-wetting k$_{r\alpha}$, Final-wetting P$_c$, K=$10^{-16}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[1][:])*dt*day/year, len(CO2[1][:])), CO2[1][:], color=[.3,0,.6], linewidth=lw, linestyle=":", label=r"Final-wetting k$_{r\alpha}$, Final-wetting P$_c$, K=$10^{-16}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[2][:])*dt*day/year, len(CO2[2][:])), CO2[2][:], color=[.8,.8,0], linewidth=lw, linestyle="-", label=r"Initial-wetting k$_{r\alpha}$, Final-wetting P$_c$, K=$7.5\times 10^{-17}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[3][:])*dt*day/year, len(CO2[3][:])), CO2[3][:], color=[.8,.8,0], linewidth=lw, linestyle=":", label=r"Final-wetting k$_{r\alpha}$, Final-wetting P$_c$, K=$7.5\times 10^{-17}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[4][:])*dt*day/year, len(CO2[4][:])), CO2[4][:], color=[1,0,.5], linewidth=lw, linestyle="-", label=r"Initial-wetting k$_{r\alpha}$, Final-wetting P$_c$, K=$5\times 10^{-17}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[5][:])*dt*day/year, len(CO2[5][:])), CO2[5][:], color=[1,0,.5], linewidth=lw, linestyle=":", label=r"Final-wetting k$_{r\alpha}$, Final-wetting P$_c$, K=$5\times 10^{-17}$ m$^2$")
plt.xlim([0,T])
plt.xlabel('t [years]')
plt.ylabel('M$_{CO_{2}}$ [kg/m$^2$]')
plt.grid()
plt.legend(loc='lower right')
matplotlib.pyplot.grid(True, which="both")
plt.savefig('fig8.eps', format='eps')
plt.show()
