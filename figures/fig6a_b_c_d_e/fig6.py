# Setting up of the system to produce fig6.
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
N = 4
ii=0
P = []
a = []
a.append(cic)
a.append(Ei)
P.append(a)
a = []
a.append(cic)
a.append(Ef)
P.append(a)
a = []
a.append(cfc)
a.append(Ei)
P.append(a)
a = []
a.append(cfc)
a.append(Ef)
P.append(a)

#Create variables used in the simulations
bashi = [0]*N
Oki = [0]*N

#Delete previous simulation files
os.system('rm -r vtk & wait')

#Set the different simulations
os.system('mkdir vtk')
for i in range(N):
    os.system('mkdir vtk/vtk-%05d & wait' % i)
    bashi[i]="%s --output-dir=vtk/vtk-%05d --wa-vtk-time-steps-file=writetimes.DATA --initial-time-step-size=.01 --max-time-step-size=86400 --max-time-step-divisions=20  --end-time=%d --enable-wa=true --beta=%E --eta=%E --ci-c=%f --cf-c=%f --ci=%f --cf=%f --ei=%f --ef=%f --lambda=%f --llambda=%f --srw=%f --srn=%f --s0w=%f --x0n=%f --k-c=%E --k=%E --phi-c=%f --phi=%f --fine-layer-bottom=%f " % (flowwa,i,T*year,1e15,0,P[i][0],P[i][0],cia,cfa,P[i][1],P[i][1],lambdaa,Lambdaa,Srw,Srn,S0w,X0n,Kc,Ka,phic,phia,h)

#Create the .bash file and run the simulations
bash=bashi[0]+" & "+bashi[1]+" & "+bashi[2]+" & "+bashi[3]+" & wait\n"
a_file = open("fig6.bash", "w")
a_file.writelines(bash)
a_file.close()
os.system("chmod u+x ./fig6.bash")
os.system('./fig6.bash')

#Check that all different simulations finished (the number of output files should be the same (+-1))
for i in range(N):
    a = []
    a_file = open("vtk/vtk-%05d/wa_ncp_ecfv.pvd" % i,"r")
    list_of_lines = a_file.readlines()
    ss=list_of_lines[-3]
    Oki[i]=int(ss[-13:-8])
print(Oki)

#Plot and save the results for fig6a_b_c_d_e.eps (the figure in the paper was produced using paraview)
sargs = dict(title_font_size=12,label_font_size=12,shadow=True,n_labels=4,italic=True,fmt="%.1f",font_family="arial",height=0.6, vertical=True, position_x=0.7, position_y=0.2,color='black')
mesh = pv.read('vtk/vtk-00000/wa_ncp_ecfv-00000.vtu')
mesh.set_active_scalars('saturation_gas')
p = pv.Plotter(shape=(1,5), window_size=(800, 250))

p.subplot(0, 0)
p.add_text("(a)", font_size=12, color='black')
p.add_mesh(mesh, scalar_bar_args=sargs, cmap="Blues", clim=[0, 0.8])
p.view_xz()
p.camera.Zoom(1)
p.show_bounds(xlabel='x [m]',zlabel='z [m]',show_xaxis=0, color='black')

p.subplot(0, 1)
p.add_text("(b)", font_size=12, color='black')
a_file = open("vtk/vtk-%05d/wa_ncp_ecfv.pvd" % 0,"r")
list_of_lines = a_file.readlines()
ss=list_of_lines[-3]
mesh = pv.read("vtk/vtk-%05d/wa_ncp_ecfv-%05d.vtu" % (0,int(ss[-13:-8])))
mesh.set_active_scalars('saturation_gas')
p.add_mesh(mesh, scalar_bar_args=sargs, cmap="Blues", clim=[0, 0.8])
p.view_xz()
p.show_bounds(xlabel='x [m]',zlabel='z [m]',show_xaxis=0, color='black')

p.subplot(0, 2)
p.add_text("(c)", font_size=12, color='black')
a_file = open("vtk/vtk-%05d/wa_ncp_ecfv.pvd" % 1,"r")
list_of_lines = a_file.readlines()
ss=list_of_lines[-3]
mesh = pv.read("vtk/vtk-%05d/wa_ncp_ecfv-%05d.vtu" % (1,int(ss[-13:-8])))
mesh.set_active_scalars('saturation_gas')
p.add_mesh(mesh, scalar_bar_args=sargs, cmap="Blues", clim=[0, 0.8])
p.view_xz()
p.show_bounds(xlabel='x [m]',zlabel='z [m]',show_xaxis=0, color='black')

p.subplot(0, 3)
p.add_text("(d)", font_size=12, color='black')
a_file = open("vtk/vtk-%05d/wa_ncp_ecfv.pvd" % 2,"r")
list_of_lines = a_file.readlines()
ss=list_of_lines[-3]
mesh = pv.read("vtk/vtk-%05d/wa_ncp_ecfv-%05d.vtu" % (2,int(ss[-13:-8])))
mesh.set_active_scalars('saturation_gas')
p.add_mesh(mesh, scalar_bar_args=sargs, cmap="Blues", clim=[0, 0.8])
p.view_xz()
p.show_bounds(xlabel='x [m]',zlabel='z [m]',show_xaxis=0, color='black')

p.subplot(0, 4)
p.add_text("(e)", font_size=12, color='black')
a_file = open("vtk/vtk-%05d/wa_ncp_ecfv.pvd" % 3,"r")
list_of_lines = a_file.readlines()
ss=list_of_lines[-3]
mesh = pv.read("vtk/vtk-%05d/wa_ncp_ecfv-%05d.vtu" % (3,int(ss[-13:-8])))
mesh.set_active_scalars('saturation_gas')
p.add_mesh(mesh, scalar_bar_args=sargs, cmap="Blues", clim=[0, 0.8])
p.view_xz()
p.show_bounds(xlabel='x [m]',zlabel='z [m]',show_xaxis=0, color='black')
p.set_background(color='white')

p.save_graphic('fig6a_b_c_d_e.eps')
p.show()
