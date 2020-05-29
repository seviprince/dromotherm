import numpy as np
import matplotlib.pyplot as plt
from dromosense.tools import sol_tridiag,Tsorties_echangeur
from dromosense.constantes import rho_eau,Cpf
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from examples.mod1D.batiment2 import Besoin,start,summerStart,summerEnd,step
#from examples.mod1D.code1D_ncouches1 import T,kelvin,qf,meteo
import cmath as math
"""
on peut importer automatiquement albedo, epsilon et sigma
from dromosense.constantes import *
"""

kelvin = 273.15
meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')
print(meteo.shape)
qf = 0.035/3600.0
T = np.loadtxt('T1d.txt')#+kelvin
print(T.shape)
input("press any key")

# Stockage thermique

# équation diff régissant la température du stockage

def F(y,t):
    """
    Le stockage a 2 entrées et 2 sorties
    Te1=Ts2 (entrée du fluide quittant le dromotherm vers le stockage)
    
    Te2=Température d'entrée de l'eau glycolée quittant la PAC pour le stockage
    
    Ts1=Tsorties_stockage(rho_eau*qf,Cpf,Te1,y,k_global)
    """
    print(t)
    i=int(t)
    Te1=Ts2_ech[i]
    Te2=Te2_sto[i]
    return (1/(m_sable*Cp_sable))*((m1*Cpf))*(Te1-Tsorties_stockage(rho_eau*qf,Cpf,Te1,y,k_global))+m2*C2*(Te2-Tsorties_stockage(m2,C2,Te2,y,k_global))
    
def Tsorties_stockage(m,Cp,Te,Tsable,k):
    """
    expliquer
    """
    return (1/(m*Cp+0.5*k))*((m*Cp-0.5*k)*Te+k*Tsable)


Ts1_ech=np.zeros(meteo.shape[0]) # Ts1_ech=Te1_sto
Ts2_ech=np.zeros(meteo.shape[0]) 
Ts1_sto=10*np.ones(meteo.shape[0])
Ts2_sto=np.zeros(meteo.shape[0])
Te2_sto=np.zeros(meteo.shape[0])
# Températures initailes de tout le système
Ts1_ech[0]=10#+kelvin
Ts2_ech[0]=10#+kelvin
Ts1_sto[0]=10#+kelvin
Ts2_sto[0]=10#+kelvin
Te2_sto[0]=10#+kelvin
# Température fluide circulant dans le dromotherm
Tsf_dro=T[:,1]

# DONNES SUR LE STOCKAGE
eff=0.8 # efficacité de l'échangeur 
q1=qf # débit du fluide circulant de l'échangeur vers le stockage
m1=rho_eau*q1
q2=0.004/60 # débit du fluide caloporteur (eau glycolée)
rho2=1040.0 # masse volumique de l'eau glycolée
m2=rho2*q2 # débit massique de l'eau glycolée
Cp2=3942.0 # capacité calorifique massique de l'eau glycolée
m_sable=70200.0 # masse de sable en kg
Cp_sable=1470.0 # capacité calorifique massqieu du sable en J/Kg.K
R1=0.0102 # Rayon intérieur du tube
R2=0.0125 # Rayon extérieur 
L_tube=5 # Longueur tube
N_tube=16  # Nombre de tube
Nu=4.36 #Nombre de tube
lambda2=15.8
lambda_tube=1.32
Rcond=np.log(R2/R1)/(2*math.pi*N_tube*L_tube*lambda_tube)
Rconv=1/(math.pi*N_tube*L_tube*lambda2)
k_global=1/(Rcond+Rconv)
m2=1000*q2
C2=4180.0 
Tsol=10.0
e_iso=0.001
S_iso=2*math.pi*2*R2*N_tube*L_tube


"""
SUMMER SIMULATION

il faut séquencer le code exactement comme celà doit se passer dans la réalité

on calcule les index permettant de boucler sur l'été
"""
i_summerStart=(summerStart-start)//step
i_summerEnd=i_summerStart+(summerEnd-summerStart)//step
print(i_summerStart)
print(i_summerEnd)
input("press any key")
for i in range (i_summerStart,i_summerEnd):

    # Températures fluides sortant de l'échangeur vers le stockage
    """
    1 = SC : source chaude = fournie par le dromotherme
    
    2 = SF : source froide = fournie par les sondes ou les corbeilles
    
    Ts1_ech: Température de sortie de SC dromothermique après passage dans l'échangeur de séparation de réseau
    
    Te2_ech: Température d'entrée de SF, du stockage vers l'échangeur de séparation de réseau = Ts1_sto ?
    
    Ts2_ech: Température de sortie de SF, de l'échangeur de séparation de réseau vers le stockage
    
    """
    Ts1_ech[i],Ts2_ech[i] = Tsorties_echangeur(Tsf_dro[i], Ts1_sto[i], rho_eau*qf, m1, Cpf, Cpf, eff)
    
    """
    
    désormais il faut calculer la température du stockage, que tu appelles Tsable
    
    comment fait-on là ? 
    
    ???
    """
    
    
plt.subplot(211)
plt.plot(Tsf_dro,label="Te1",color="red")
plt.plot(Ts1_ech,label="Ts1",color="purple")
plt.legend()
plt.subplot(212)
plt.plot(Ts1_sto,label="Te2",color="blue")
plt.plot(Ts2_ech,label="Ts2",color="orange")
plt.legend()
plt.show()
    

#Température du stockage
Tsable = odeint(F,10,meteo[0:-10,0])



"""
# Pompe à Chaleur

COP=3
for i in range (0,meteo.shape[0]-1):
    
    Te2_sto[i]=Ts2_sto[i]-COP*Besoin[i]/((COP-1)*m2*Cp2)
    Ts1_sto[i]=Tsorties_stockage(m1,Cpf,Ts2_ech[i],Tsable[i],k_global)
    Ts2_sto[i]=Tsorties_stockage(m2,Cp2,Te2_sto[i],Tsable[i],k_global)
"""
    









