# Example to run wa.
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
import pyvista as pv

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
tch=1e7             #Characterictic time [s]
b1=1e7              #Dynamic parameter (capillary pressure) [-]
b2=1.8              #Dynamic parameter (capillary pressure) [-]
n1=4.999e3          #Dynamic parameter (relative permeability) [-]
n2=5e-1             #Dynamic parameter (relative permeability) [-]
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
T=100.              #Total simulation time [y]
dt=365.             #Time step to print the results [d]
C=1e-5              #Pore-scale parameter [-]

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

#Create variables used in the simulations
CO2 = []

#Delete previous simulation files
os.system('rm -r vtk & wait')

#Set the simulation
os.system('mkdir vtk')
bash="%s --output-dir=vtk --wa-vtk-time-steps-file=writetimes.DATA --initial-time-step-size=.01 --max-time-step-size=2592000 --max-time-step-divisions=20  --end-time=%d --enable-wa=true --tch=%E --beta=%E --eta=%E --ci-c=%f --cf-c=%f --ci=%f --cf=%f --ei=%f --ef=%f --lambda=%f --llambda=%f --srw=%f --srn=%f --s0w=%f --x0n=%f --k-c=%E --k=%E --phi-c=%f --phi=%f --fine-layer-bottom=%f & wait\n" % (flowwa,T*year,tch,b1*C**b2,-n1*C+n2,cic,cfc,cia,cfa,Ei,Ef,lambdaa,Lambdaa,Srw,Srn,S0w,X0n,Kc,Ka,phic,phia,h)

#Create the .bash file and run the simulation
a_file = open("wa.bash", "w")
a_file.writelines(bash)
a_file.close()
os.system("chmod u+x ./wa.bash")
os.system('./wa.bash')

#Compute the CO2 in the caprock from the vtk output files
mesh = meshio.read("vtk/wa_ncp_ecfv-00000.vtu")
for row in mesh.cell_data['saturation_gas']:
    sn=row
inx=np.where(sn<1e-10)
a = []
a_file = open("vtk/wa_ncp_ecfv.pvd","r")
list_of_lines = a_file.readlines()
ss=list_of_lines[-3]
for j in range(int(ss[-13:-8])+1):
    mesh = meshio.read("vtk/wa_ncp_ecfv-%05d.vtu" % j)
    for row in mesh.cell_data['saturation_gas']:
        snd=row
    a_file.close()
    a.append(sum(rhon*phic*snd[inx]))
CO2.append(a)

#Plot the non-wetting saturation at the end of the simulation and save it as snfinal.eps
sargs = dict(title_font_size=12,label_font_size=12,shadow=True,n_labels=4,italic=True,fmt="%.1f",font_family="arial",height=0.6, vertical=True, position_x=0.7, position_y=0.2,color='black')
p = pv.Plotter(shape=(1, 1), window_size=(200, 250))
a_file = open("vtk/wa_ncp_ecfv.pvd","r")
list_of_lines = a_file.readlines()
ss=list_of_lines[-3]
mesh = pv.read("vtk/wa_ncp_ecfv-%05d.vtu" % int(ss[-13:-8]))
mesh.set_active_scalars('saturation_gas')
p.add_mesh(mesh, scalar_bar_args=sargs, cmap="Blues", clim=[0, 1-Srw])
p.view_xz()
p.camera.Zoom(1)
p.show_bounds(xlabel='x [m]',zlabel='z [m]',show_xaxis=0, color='black')
p.set_background(color='white')
p.save_graphic('snfinal.eps')
p.show()

#Plot the CO2 mass/area in the caprock and save it as co2final.eps
lw=5
plt.figure(figsize=(5, 5), dpi=80)
plt.rc('font', size=9)
plt.rc('legend', fontsize=9)
axes=plt.subplot(1, 1, 1)
plt.plot(np.linspace(0, len(CO2[0][:])*dt*day/year, len(CO2[0][:])), CO2[0][:], color=[1,.5,0], linewidth=lw, linestyle=":")
plt.xlim([0,T])
plt.xlabel('t [years]')
plt.ylabel('M$_{CO_{2}}$ [kg/m$^2$]')
plt.grid()
matplotlib.pyplot.grid(True, which="both")
plt.legend(loc='upper left')
plt.savefig('co2final.eps', format='eps')
plt.show()
