from __future__ import division
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.pyplot import *
from scipy.stats import norm
from scipy.interpolate import interp1d
import math as m

#AJHC
#Umean=58.74
#Uco=3.2
#Uair=3.3
############################################## INPUT DATA #################################################################################
Re=6191890.08267569           		             		# Jet Reynolds number based on the nozzle diameter
R=0.0098*0.5       			     			# Radius of the inlet
h_otw=R
Umean=6.22
x1 = np.loadtxt('inlet_coords.txt', unpack=True)        	# The cell centers coordinates
x2 = np.loadtxt('inlet_coords.txt', unpack=True)-(h_otw) 	# The radial coordinates
print(x2)
case_name='wycinek_6.22'
#######################################################################################################################################
################################################### Turbulent velocity profile ########################################################
############################### Maximum velocity
n=(2.1*m.log(Re,10)-1.9)
turb_coeff=2*n**2/((n+1)*(2*n+1))
print(turb_coeff)
print(n)
Umax=Umean/turb_coeff
print( 'Umax:')
print (Umax)
############################### Velocty profile definition - "power law"
def Ujet(r):
    if r<0:
        r=-r
    return Umax*(1-(r/R))**(1/n)
Xj=[]
for i in np.linspace(-R, R, 2001):
    Xj.append(i)
Ujet_turb = []
for i in range(len(Xj)):
    Ujet_turb.append(Ujet(Xj[i]))
Ujet_turb = np.array(Ujet_turb)

Ujet_chart = []
for i in range(len(x2)):
    Ujet_chart.append(Ujet(x2[i]))
Ujet_chart = np.array(Ujet_chart)

############################### Get the cell centers coordinates values
##x1 = np.loadtxt('ccy_jet', unpack=True)
########################################## plot and write to file #####################################################################
figure()
ylabel(r'$\widetilde{U}$ $[\frac{m}{s}]$', fontsize='large')
xlabel('r [m]', fontsize='large')
tick_params(axis='x', labelsize=14)
tick_params(axis='y', labelsize=14)
plot(x1, Ujet_chart, '--', color='blue', label='turbulent U profile')
plot(x1, Ujet_chart, 'o', color='blue', label='Cells values')
plot([-R+h_otw,-R+h_otw], [0,1.2*Umax], '-', color='black', label='duct_wall_1')
plot([R+h_otw,R+h_otw], [0,1.2*Umax], '-', color='black', label='duct_wall_2')
legend(loc=3, prop={'size':14})
grid()
ylim(0, 1.2*Umax)
savefig(str(case_name)+'_inlet_U.pdf')
#show()

file = open(str(case_name)+"_U_inlet.txt", "w")
for i in range(len(x1)):
    file.write("(" + str(Ujet_chart[i]) + " " + "0" + " " + "0)" + "\n")
file.write("\n")
for i in range(len(x1)):
    file.write(str(Ujet_chart[i]) + "\n")
file.close()
#######################################################################################################################################
##################################################### du/dr ###########################################################################
dr=np.diff(Xj)[0]
dudr = abs(np.diff(Ujet_turb)/dr)
del Xj[-1]
#function du(r)/dr
F_dudr = interp1d(Xj, dudr, kind='cubic')

#######################################################################################################################################
################################################### Turbulent kinetic energy profile ##################################################
#k list and k(r) function
k=[]
for i in range(0, len(Xj)):
    k.append(0.025*2*R*Re**(-0.125)*Ujet(Xj[i])*F_dudr(Xj[i]))
F_k = interp1d(Xj, k, kind='cubic')

########################################## plot and write to file #####################################################################
figure()
ylabel(r'$k$ $[\frac{m^2}{s^2}]$', fontsize='large')
xlabel('r [m]', fontsize='large')
tick_params(axis='x', labelsize=14)
tick_params(axis='y', labelsize=14)
ylim(0, 1.2*F_k(max(x2)))
#plot(Xj, k, '--', color='blue', label='k profile')
plot(x1, F_k(x2), '--', color='blue', label='turbulent k profile')
plot(x1, F_k(x2), 'o', color='blue', label='Cells values')
plot([-R+h_otw,-R+h_otw], [0,1.2*F_k(max(x2))], '-', color='black', label='duct_wall_1')
plot([R+h_otw,R+h_otw], [0,1.2*F_k(max(x2))], '-', color='black', label='duct_wall_2')
legend(loc=2, prop={'size':14})
grid()
savefig(str(case_name)+'_inlet_k.pdf')
#show()

file = open(str(case_name)+"_k_inlet.txt", "w")
for i in range(len(x2)):
    file.write(str(F_k(x2[i])) + "\n")
file.write("\n")
for i in range(len(x2)):
    file.write(str(F_k(x2[i])) + "\n ")
file.close()

#######################################################################################################################################
############################################ Turbulent kinetic energy dissipation rate profile ########################################
##Epsilon list and eps(r) function
eps=[]
for i in range(0, len(Xj)):
    eps.append(0.0075*2*R*Re**(-0.125)*Ujet(Xj[i])*F_dudr(Xj[i])**2)
F_eps = interp1d(Xj, eps, kind='cubic')

########################################## plot and write to file #####################################################################
figure()
ylabel(r'$\epsilon$ [$\frac{m^2}{s^3}$]', fontsize='large')
xlabel('r [m]', fontsize='large')
tick_params(axis='x', labelsize=14)
tick_params(axis='y', labelsize=14)
ylim(0, 1.2*F_eps(max(x2)))
plot(x1, F_eps(x2), '--', color='blue', label=r'$\epsilon$ profile')
plot(x1, F_eps(x2), 'o', color='blue', label='Cells values')
plot([-R+h_otw,-R+h_otw], [0,1.2*F_eps(max(x2))], '-', color='black', label='duct_wall_1')
plot([R+h_otw,R+h_otw], [0,1.2*F_eps(max(x2))], '-', color='black', label='duct_wall_2')
legend(loc=2, prop={'size':14})
grid()
savefig(str(case_name)+'_inlet_eps.pdf')
#show()

file = open(str(case_name)+"_eps_inlet.txt", "w")
for i in range(len(x1)):
    file.write(str(F_eps(x2[i])) + "\n")
file.write("\n")
for i in range(len(x1)):
    file.write(str(F_eps(x2[i])) + "\n ")
file.close()
#######################################################################################################################################
############################################ Omega profile ############################################################################
##Epsilon list and eps(r) function
omega=[]
for i in range(0, len(x2)):
    omega.append((F_eps(x2[i])/(0.09*F_k(x2[i]))))
    #print(omega[i])

########################################## plot and write to file #####################################################################
figure()
ylabel(r'$\omega$ [$\frac{m^2}{s^3}$]', fontsize='large')
xlabel('r [m]', fontsize='large')
tick_params(axis='x', labelsize=14)
tick_params(axis='y', labelsize=14)
ylim(-500, 1.2*max(omega))
plot(x1, omega, '--', color='blue', label=r'$\omega$ profile')
plot(x1, omega, 'o', color='blue', label='Cells values')
plot([-R+h_otw,-R+h_otw], [0,1.2*max(omega)], '-', color='black', label='duct_wall_1')
plot([R+h_otw,R+h_otw], [0,1.2*max(omega)], '-', color='black', label='duct_wall_2')
legend(loc=2, prop={'size':14})
grid()
savefig(str(case_name)+'_inlet_omega.pdf')
#show()

file = open(str(case_name)+"_omega_inlet.txt", "w")
for i in range(len(x1)):
    file.write(str(omega[i]) + "\n")
file.write("\n")
for i in range(len(x1)):
    file.write(str(omega[i]) + "\n ")
file.close()
############################################# END JET INLET BOUNDARY CONDITIONS #########################################################




