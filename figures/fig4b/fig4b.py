# Setting up of the system to produce fig4b.
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
dt=73.              #Time step to print the results [d]

#Set the values for the different simulations
N=5
M=3
C=[1e-5,2.5e-5,5e-5,7.5e-5,1e-4]
q=[1e-7,5e-8,2.5e-8]
K=[1e-10,5e-11,2.5e-11]
phi=[0.1,0.2,0.4]

#Define conversion variables
milli=1e-3          #[-]
darcy=9.8692e-13    #[m^2]
cp=1e-3             #[Pa s]
day=86400           #[s]
atm=101325          #[Pa]

#Define the cases for the numerical studies
ii=0
H = []
for j in range(N):
    for i in range(M):
        a = []
        a.append(q[i])
        a.append(K[0])
        a.append(phi[0])
        a.append(C[j])
        H.append(a)
        ii=ii+1
    for i in range(M):
        a = []
        a.append(q[0])
        a.append(K[i])
        a.append(phi[0])
        a.append(C[j])
        H.append(a)
        ii=ii+1
    for i in range(M):
        a = []
        a.append(q[0])
        a.append(K[0])
        a.append(phi[i])
        a.append(C[j])
        H.append(a)
        ii=ii+1

#Create variables used in the simulations
bash = ['1']*(ii*2)
bashi = ['1']*ii
bashd = ['1']*ii
Oki = [0]*ii
Okd = [0]*ii
Ca = [0]*ii
sfld = [0]*ii
xin = [0]*ii
xdyn = [0]*ii
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
list_of_lines[163] = "%d*%f /\n" % (int(T/dt),dt)

#Set the different simulations
os.system('mkdir vtk & wait')
os.system('mkdir WACASES & wait')
for i in range(ii):
    list_of_lines[61] = "%d*%f /\n" % (nx,H[i][2])
    list_of_lines[64] = "%d*%f /\n" % (nx,(H[i][1]/(milli*darcy)))
    list_of_lines[67] = "%d*%f /\n" % (nx,(H[i][1]/(milli*darcy)))
    list_of_lines[70] = "%d*%f /\n" % (nx,(H[i][1]/(milli*darcy)))
    list_of_lines[147] = "'INJE01' 'WATER' 'OPEN' 'RATE' %E  1* 1E8 /\n" % (H[i][0]*day)
    list_of_lines[157] = "'INJE01' 'OIL' 'OPEN' 'RATE' %E  1* 1E8 /\n" % (H[i][0]*day)
    a_file = open("WACASES/WACASEI-%05d.DATA" % i, "w")
    a_file.writelines(list_of_lines)
    a_file.close()
    a_file = open("WACASES/WACASED-%05d.DATA" % i, "w")
    a_file.writelines(list_of_lines)
    a_file.close()
    os.system('mkdir vtk/vtk-%05d & wait' % i)
    bashi[i]="%s WACASES/WACASEI-%05d.DATA --output-dir=vtk/vtk-%05d --enable-vtk-output=true --enable-ecl-output=false --initial-time-step-in-days=0.0001 --solver-max-restarts=20 --solver-max-time-step-in-days=.1 --enable-wa=true --beta=%E --eta=%E --ei=%f --ef=%f --ci=%E --cf=%E --lambda=%f --llambda=%f --srw=%f" % (flowpath,i,i,1E15,0,Ei,Ei,ci,ci,lambdaa,Lambdaa,Srw)
    bashd[i]="%s WACASES/WACASED-%05d.DATA --output-dir=vtk/vtk-%05d --enable-vtk-output=true --enable-ecl-output=false --initial-time-step-in-days=0.0001 --solver-max-restarts=20 --solver-max-time-step-in-days=.1 --enable-wa=true --beta=%E --eta=%E --ei=%f --ef=%f --ci=%E --cf=%E --lambda=%f --llambda=%f --srw=%f" % (flowpath,i,i,1E15,-n1*H[i][3]+n2,Ei,Ef,ci,ci,lambdaa,Lambdaa,Srw)

#Create the .bash file and run the simulations
j=0
for i in range(int(ii/5)):
    bash[j]=bashi[5*i]+" & "+bashi[5*i+1]+" & "+bashi[5*i+2]+" & "+bashi[5*i+3]+" & "+bashi[5*i+4]+" & wait\n"
    j=j+1
for i in range(int(ii/5)):
    bash[j]=bashd[5*i]+" & "+bashd[5*i+1]+" & "+bashd[5*i+2]+" & "+bashd[5*i+3]+" & "+bashd[5*i+4]+" & wait\n"
    j=j+1
a_file = open("fig4b.bash", "w")
a_file.writelines(bash)
a_file.close()
os.system("chmod u+x ./fig4b.bash")
os.system('./fig4b.bash')

#Obtain the variables from the vtk output files
for i in range(ii):
    a_file = open("vtk/vtk-%05d/WACASEI-%05d.pvd" % (i,i),"r")
    list_of_lines = a_file.readlines()
    ss=list_of_lines[-3]
    Oki[i]=int(ss[-13:-8])
    mesh = meshio.read("vtk/vtk-%05d/WACASEI-%05d-%05d.vtu" % (i,i,int(ss[-13:-8])))
    for row in mesh.cell_data['saturation_oil']:
        sni=row
    a_file.close()
    a_file = open("vtk/vtk-%05d/WACASED-%05d.pvd" % (i,i),"r")
    list_of_lines = a_file.readlines()
    ss=list_of_lines[-3]
    Okd[i]=int(ss[-13:-8])
    mesh = meshio.read("vtk/vtk-%05d/WACASED-%05d-%05d.vtu" % (i,i,int(ss[-13:-8])))
    for row in mesh.cell_data['saturation_oil']:
        snd=row
    a_file.close()
    Ca[i]= (H[i][0]*mun)/(H[i][1]*H[i][2]*ci)
    xin[i]=(nx-len(sni[tuple([sni<1e-12])])-.5)*dx
    xdyn[i]=(nx-len(snd[tuple([snd<1e-12])])-.5)*dx
    sfld[i]=(xin[i]-xdyn[i])/xin[i]

#Check that all different simulations finished (the number of output files should be the same)
print(Oki)
print(Okd)

#Plot the results and save them as fig4b.eps
lw=1
plt.figure(figsize=(4.5, 4), dpi=512)
plt.rc('font', size=9)
axes=plt.subplot(1, 1, 1)
axes.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.plot(Ca[0:3], sfld[0:3], color=[0,0,0], linewidth=lw, linestyle="--", marker="v", label="q (C=10$^{-5}$)")
plt.plot(Ca[3:6], sfld[3:6], color=[1,0,0], linewidth=lw, linestyle="--", marker="v", label="K (C=10$^{-5}$)")
plt.plot(Ca[6:9], sfld[6:9], color=[0,0,1], linewidth=lw, linestyle="--", marker="v", label="$\phi$ (C=10$^{-5}$)")
plt.plot(Ca[9:12], sfld[9:12], color=[0,0,0], linewidth=lw, linestyle="--", marker="o", label=r"q (C=2.5 $\times 10^{-5}$)")
plt.plot(Ca[12:15], sfld[12:15], color=[1,0,0], linewidth=lw, linestyle="--", marker="o", label=r"K (C=2.5 $\times 10^{-5}$)")
plt.plot(Ca[15:18], sfld[15:18], color=[0,0,1], linewidth=lw, linestyle="--", marker="o", label=r"$\phi$ (C=2.5 $\times 10^{-5}$)")
plt.plot(Ca[18:21], sfld[18:21], color=[0,0,0], linewidth=lw, linestyle="--", marker="x", label=r"q (C=5 $\times 10^{-5}$)")
plt.plot(Ca[21:24], sfld[21:24], color=[1,0,0], linewidth=lw, linestyle="--", marker="x", label=r"K (C=5 $\times 10^{-5}$)")
plt.plot(Ca[24:27], sfld[24:27], color=[0,0,1], linewidth=lw, linestyle="--", marker="x", label=r"$\phi$ (C=5 $\times 10^{-5}$)")
plt.plot(Ca[27:30], sfld[27:30], color=[0,0,0], linewidth=lw, linestyle="--", marker="s", label=r"q (C=7.5 $\times 10^{-5}$)")
plt.plot(Ca[30:33], sfld[30:33], color=[1,0,0], linewidth=lw, linestyle="--", marker="s", label=r"K (C=7.5 $\times 10^{-5}$)")
plt.plot(Ca[33:36], sfld[33:36], color=[0,0,1], linewidth=lw, linestyle="--", marker="s", label=r"$\phi$ (C=7.5 $\times 10^{-5}$)")
plt.plot(Ca[36:39], sfld[36:39], color=[0,0,0], linewidth=lw, linestyle="--", marker="^", label="q (C=10$^{-4}$)")
plt.plot(Ca[39:42], sfld[39:42], color=[1,0,0], linewidth=lw, linestyle="--", marker="^", label="K (C=10$^{-4}$)")
plt.plot(Ca[42:45], sfld[42:45], color=[0,0,1], linewidth=lw, linestyle="--", marker="^", label="$\phi$ (C=10$^{-4}$)")
plt.xlim([1e-6,3e-4])
plt.ylim([-.001,0.05])
plt.xscale('log')
plt.xlabel('Ca [-]')
plt.ylabel('SFLD [-]')
plt.grid()
plt.title('(b)')
matplotlib.pyplot.grid(True, which="both")
plt.legend(loc='best', prop={'size': 8})
plt.savefig('fig4b.eps', format='eps')
plt.show()
