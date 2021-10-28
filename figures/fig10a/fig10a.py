# Setting up of the system to produce fig10a.
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
T=100               #Total simulation time [y]
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

#Define the cases for the simulations
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
N = len(P) + 1

#Create variables used in the simulations
bash = ['1']*N
bashd = ['1']*N
Okd = [0]*N
CO2 = []

#Delete previous simulation files
os.system('rm -r vtk & wait')

#Set the different simulations
os.system('mkdir vtk')
for i in range(N):
    os.system('mkdir vtk/vtk-%05d & wait' % i)
    if (i==0):
        bashd[i]="%s --output-dir=vtk/vtk-%05d --wa-vtk-time-steps-file=writetimes.DATA --initial-time-step-size=.01 --max-time-step-size=2592000 --max-time-step-divisions=20  --end-time=%d --enable-wa=true --tch=%E --beta=%E --eta=%E --ci-c=%f --cf-c=%f --ci=%f --cf=%f --ei=%f --ef=%f --lambda=%f --llambda=%f --srw=%f --srn=%f --s0w=%f --x0n=%f --k-c=%E --k=%E --phi-c=%f --phi=%f --fine-layer-bottom=%f " % (flowwa,i,T*year,tch,1e15,0,cfc,cfc,cia,cfa,Ei,Ei,lambdaa,Lambdaa,Srw,Srn,S0w,X0n,Kc,Ka,phic,phia,h)
    else:
        bashd[i]="%s --output-dir=vtk/vtk-%05d --wa-vtk-time-steps-file=writetimes.DATA --initial-time-step-size=.01 --max-time-step-size=2592000 --max-time-step-divisions=20  --end-time=%d --enable-wa=true --tch=%E --beta=%E --eta=%E --ci-c=%f --cf-c=%f --ci=%f --cf=%f --ei=%f --ef=%f --lambda=%f --llambda=%f --srw=%f --srn=%f --s0w=%f --x0n=%f --k-c=%E --k=%E --phi-c=%f --phi=%f --fine-layer-bottom=%f " % (flowwa,i,T*year,tch,b1*P[i-1][2]**b2,0,P[i-1][0],cfc,cia,cfa,Ei,Ei,lambdaa,Lambdaa,Srw,Srn,S0w,X0n,P[i-1][1],Ka,phic,phia,h)
bashd[11]="%s --output-dir=vtk/vtk-%05d --wa-vtk-time-steps-file=writetimes.DATA --initial-time-step-size=.01 --max-time-step-size=86400 --max-time-step-divisions=20  --end-time=%d --enable-wa=true --tch=%E --beta=%E --eta=%E --ci-c=%f --cf-c=%f --ci=%f --cf=%f --ei=%f --ef=%f --lambda=%f --llambda=%f --srw=%f --srn=%f --s0w=%f --x0n=%f --k-c=%E --k=%E --phi-c=%f --phi=%f --fine-layer-bottom=%f " % (flowwa,11,T*year,tch,b1*P[10][2]**b2,0,P[10][0],cfc,cia,cfa,Ei,Ei,lambdaa,Lambdaa,Srw,Srn,S0w,X0n,P[10][1],Ka,phic,phia,h)
bashd[12]="%s --output-dir=vtk/vtk-%05d --wa-vtk-time-steps-file=writetimes.DATA --initial-time-step-size=.01 --max-time-step-size=86400 --max-time-step-divisions=20  --end-time=%d --enable-wa=true --tch=%E --beta=%E --eta=%E --ci-c=%f --cf-c=%f --ci=%f --cf=%f --ei=%f --ef=%f --lambda=%f --llambda=%f --srw=%f --srn=%f --s0w=%f --x0n=%f --k-c=%E --k=%E --phi-c=%f --phi=%f --fine-layer-bottom=%f " % (flowwa,12,T*year,tch,b1*P[11][2]**b2,0,P[11][0],cfc,cia,cfa,Ei,Ei,lambdaa,Lambdaa,Srw,Srn,S0w,X0n,P[11][1],Ka,phic,phia,h)
bashd[13]="%s --output-dir=vtk/vtk-%05d --wa-vtk-time-steps-file=writetimes.DATA --initial-time-step-size=.01 --max-time-step-size=100000 --max-time-step-divisions=20  --end-time=%d --enable-wa=true --tch=%E --beta=%E --eta=%E --ci-c=%f --cf-c=%f --ci=%f --cf=%f --ei=%f --ef=%f --lambda=%f --llambda=%f --srw=%f --srn=%f --s0w=%f --x0n=%f --k-c=%E --k=%E --phi-c=%f --phi=%f --fine-layer-bottom=%f " % (flowwa,13,T*year,tch,b1*P[12][2]**b2,0,P[12][0],cfc,cia,cfa,Ei,Ei,lambdaa,Lambdaa,Srw,Srn,S0w,X0n,P[12][1],Ka,phic,phia,h)

#Create the .bash file and run the simulations
j=0
for i in range(int(N/5)):
    bash[j]=bashd[5*i]+" & "+bashd[5*i+1]+" & "+bashd[5*i+2]+" & "+bashd[5*i+3]+" & "+bashd[5*i+4]+" & wait\n"
    j=j+1
bash[j]=bashd[-1]+" & wait\n"
a_file = open("fig10a.bash", "w")
a_file.writelines(bash)
a_file.close()
os.system("chmod u+x ./fig10a.bash")
os.system('./fig10a.bash')

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
    Okd[i]=int(ss[-13:-8])
    for j in range(int(ss[-13:-8])+1):
        mesh = meshio.read("vtk/vtk-%05d/wa_ncp_ecfv-%05d.vtu" % (i,j))
        for row in mesh.cell_data['saturation_gas']:
            snd=row
        a_file.close()
        a.append(sum(rhon*phic*snd[inx]))
    CO2.append(a)

#Check that all different simulations finished (the number of output files should be the same (+-1))
print(Okd)

#Plot and save the results for fig10a.eps
lw=3
nn=3
plt.figure(figsize=(8, 6), dpi=80)
plt.rc('font', size=12)
plt.rc('legend', fontsize=7)

axes=plt.subplot(1, 1, 1)
plt.plot(np.linspace(0, len(CO2[0][:])*dt*day/year, len(CO2[0][:])), CO2[0][:], linewidth=lw, linestyle=":", color=[0.3,0,0.6], label=r"K=$10^{-16}$ m$^2$, Initial-wetting k$_{r\alpha}$, Final-wetting P$_c$")
plt.plot(np.linspace(0, len(CO2[1][:])*dt*day/year, len(CO2[1][:])), CO2[1][:], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="o",markevery=nn, label=r"K=$10^{-16}$ m$^2$, c$^i$=$10^4$ Pa, C=$10^{-5}$")
plt.plot(np.linspace(0, len(CO2[2][:])*dt*day/year, len(CO2[2][:])), CO2[2][:], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="s",markevery=nn, label=r"K=$7.5\times 10^{-17}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[3][:])*dt*day/year, len(CO2[3][:])), CO2[3][:], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="d",markevery=nn, label=r"K=$5\times 10^{-17}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[4][:])*dt*day/year, len(CO2[4][:])), CO2[4][:], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="P",markevery=nn, label=r"$c^i=5\times 10^4$ Pa")
plt.plot(np.linspace(0, len(CO2[5][:])*dt*day/year, len(CO2[5][:])), CO2[5][:], mec=[0,0,0], mfc=[0,0,0], linestyle="", marker="v",markevery=nn, label=r"$c^i=10^5$ Pa")
plt.plot(np.linspace(0, len(CO2[6][:])*dt*day/year, len(CO2[6][:])), CO2[6][:], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="o",markevery=nn, label=r"K=$10^{-16}$ m$^2$, c$^i$=$10^4$ Pa, C=$2.5\times 10^{-5}$")
plt.plot(np.linspace(0, len(CO2[7][:])*dt*day/year, len(CO2[7][:])), CO2[7][:], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="s",markevery=nn, label=r"K=$7.5\times 10^{-17}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[8][:])*dt*day/year, len(CO2[8][:])), CO2[8][:], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="d",markevery=nn, label=r"K=$5\times 10^{-17}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[9][:])*dt*day/year, len(CO2[9][:])), CO2[9][:], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="P",markevery=nn, label=r"$c^i=5\times 10^4$ Pa")
plt.plot(np.linspace(0, len(CO2[10][:])*dt*day/year, len(CO2[10][:])), CO2[10][:], mec=[1,0,0], mfc=[1,0,0], linestyle="", marker="v",markevery=nn, label=r"$c^i=10^5$ Pa")
plt.plot(np.linspace(0, len(CO2[11][:])*dt*day/year, len(CO2[11][:])), CO2[11][:], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="o",markevery=nn, label=r"K=$10^{-16}$ m$^2$, c$^i$=$10^4$ Pa, C=$5\times 10^{-5}$")
plt.plot(np.linspace(0, len(CO2[12][:])*dt*day/year, len(CO2[12][:])), CO2[12][:], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="s",markevery=nn, label=r"K=$7.5\times 10^{-17}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[13][:])*dt*day/year, len(CO2[13][:])), CO2[13][:], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="d",markevery=nn, label=r"K=$5\times 10^{-17}$ m$^2$")
plt.plot(np.linspace(0, len(CO2[14][:])*dt*day/year, len(CO2[14][:])), CO2[14][:], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="P",markevery=nn, label=r"$c^i=5\times 10^4$ Pa")
plt.plot(np.linspace(0, len(CO2[15][:])*dt*day/year, len(CO2[15][:])), CO2[15][:], mec=[0,0,1], mfc=[0,0,1], linestyle="", marker="v",markevery=nn, label=r"$c^i=10^5$ Pa")

plt.xlim([0,T])
plt.ylim([0,250])
plt.xlabel('t [years]')
plt.ylabel('M$_{CO_{2}}$ [kg/m$^2$]')
plt.grid()
matplotlib.pyplot.grid(True, which="both")
plt.legend(loc='upper left')
plt.savefig('fig10a.eps', format='eps')
plt.show()

#Save the data for fig11a
os.system('mkdir data-for-fig11a')
for i in range(N):
    np.savetxt('data-for-fig11a/wa-fig10a-co2-%05d.csv' % i, CO2[i][:], delimiter=',')
    np.savetxt('data-for-fig11a/wa-fig10a-time-%05d.csv' % i, np.linspace(0, len(CO2[i][:])*dt*day/year, len(CO2[i][:])), delimiter=',')
