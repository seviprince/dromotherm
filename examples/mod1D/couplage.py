import numpy as np
import matplotlib.pyplot as plt
from dromosense.tools import *
from dromosense.constantes import rho_eau,Cpf
from scipy.integrate import odeint
#from scipy.integrate import solve_ivp
import cmath as math

verbose = False

start = 1483232400
summerStart = 1496278800
#summerEnd = 1506819600 : On s'était trompé sur la période de récupération de chaleur que Monsieur Frédéric a envoyée..c'est du 1er juin au 31 août et non au 30 septembre
summerEnd=1504141200
step=3600

# pour l'instant tous les calculs sont en degrés
kelvin = 273.15

meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')
print(meteo.shape)
T = np.loadtxt('T1d.txt')#+kelvin
print(T.shape)
input("press any key")


def F(y,t):
    """
    y : Tsable = Tstockage
    
    t : time index
    
    result : dy/dt = dTsable/dt = dTstockage/dt
    """
    i = int(t)
    if verbose:
        print("we have t={} and y={}".format(i,y))
    
    #Tsor_sto[i] = ( k * y / (msto * cpsto + k/2) - coeff * eff * Tsor_dro[i] ) / a
    
    Tsor_sto[i] = ( k * y - B * Tsor_dro[i] ) / A
    
    Tinj_sto[i] = Tsor_sto[i] + coeff * eff * (Tsor_dro[i] - Tsor_sto[i])
    
    Tinj_dro[i] = Tsor_dro[i] - eff * (Tsor_dro[i] - Tsor_sto[i])
    
    der=msto * cpsto * (Tinj_sto[i] - Tsor_sto[i]) / (m_sable * Cp_sable)
    
    if verbose:
        print("dTsable/dt is {}".format(der))
    
    return der


# température d'entrée et de sortie du fluide dans le stockage
Tinj_sto=np.zeros(meteo.shape[0])
# température de sortie du fluide après transit dans le stockage
Tsor_sto=np.zeros(meteo.shape[0])

# température d'injection du fluide dans le dromotherme
Tinj_dro=10*np.ones(meteo.shape[0])
# température de sortie de fluide après collecte de chaleur dans le dromotherme
Tsor_dro=T[:,1]

"""
# DONNES SUR LE STOCKAGE
eff=0.8 # efficacité de l'échangeur 
q1=qf # débit du fluide circulant de l'échangeur vers le stockage
m1=rho_eau*q1
q2=0.004/60 # débit du fluide caloporteur (eau glycolée)
rho2=1040.0 # masse volumique de l'eau glycolée
m2=rho2*q2 # débit massique de l'eau glycolée
Cp2=3942.0 # capacité calorifique massique de l'eau glycolée
"""

"""
massif de stockage
"""
#m_sable=70200.0 # masse de sable en kg
m_sable=100
Cp_sable=1470.0 # capacité calorifique massique du sable en J/Kg.K

"""
paramètres définissant le système géothermique qui équipant le stockage
"""
R1=0.0102 # Rayon intérieur du tube
R2=0.0125 # Rayon extérieur 
L_tube=5 # Longueur tube en m
N_tube=16  # Nombre de tubes
# conductivités en W/(K.m)
lambda2=15.8
lambda_tube=1.32
Nu=4.36 # Nombre de Nusselt de l'écoulement du fluide à l'intérieur des tubes
Rcond=np.log(R2/R1)/(2*math.pi*N_tube*L_tube*lambda_tube)
Rconv=1/(math.pi*Nu*L_tube*lambda2)
# k exprimé en W/K
k_global=1/(Rcond+Rconv)

print("le k du système géothermique vaut {} W/K".format(k_global))


# efficacité de l'échangeur
eff = 0.8
# débits dans chacune des parties de l'échangeur
qdro = 0.035/3600 # m3/s
qsto = 1.2*qdro
# rho_eau en provenance du fichier constantes est expérimée en kg/m3
# Cpf en provenance du fichier constantes est exprimée en J/(kg.K)
mdro = qdro * rho_eau
msto = qsto * rho_eau
cpdro = Cpf
cpsto = Cpf

k = k_global
    
coeff = (mdro * cpdro) / (msto * cpsto)

#print("(msto cpsto -k/2) / (msto cpsto + k/2) est égal à {}".format((msto * cpsto - k/2) / (msto * cpsto + k/2))) 
    
#a = 1 - coeff * eff - (msto * cpsto - k/2) / (msto * cpsto + k/2)

#print("coeff vaut {} et a vaut {}".format(coeff, a))

B = (msto * cpsto +k/2) * coeff * eff

A = k - B

print("coeff vaut {} B vaut {} W/K et A vaut {} W/K".format(coeff, B, A))

"""
m2=1000*q2
C2=4180.0 
Tsol=10.0
e_iso=0.001
S_iso=2*math.pi*2*R2*N_tube*L_tube
"""

"""
SUMMER SIMULATION

il faut séquencer le code exactement comme celà doit se passer dans la réalité

on calcule les index permettant de boucler sur l'été
"""
i_summerStart=(summerStart-start)//step
i_summerEnd=i_summerStart+(summerEnd-summerStart)//step
print("nous allons simuler la récolte énergétique entre les heures {} et {}".format(i_summerStart,i_summerEnd))
input("press any key")
    
#Température du stockage/sable
Tsable = odeint(F,10,meteo[i_summerStart:i_summerEnd,0])

"""
WINTER SIMULATION 

agenda de présence et besoin en chauffage
"""
nbpts=meteo.shape[0]
schedule=np.array([[8,17],[8,17],[8,17],[8,17],[8,17],[-1,-1],[-1,-1]])
agenda=basicAgenda(nbpts,step,start,summerStart,summerEnd,schedule=schedule)

Tconsigne=19
Rm=8.24E-02 # Résistance thermique des murs (K/W)
Ri=1.43E-03 # Résistance superficielle intérieure
Rf=0.034 # Résistance due aux infiltrations+vitre et au renouvellement d'air 

Sm=49.5
FSm=0.0048
Sv=3.5
FSv=0.8
"""
hypothèse : la maison fait 5 mètres de large sur 5 mètres de long
vu du ciel c'est donc un carré de 25m2 de surface
"""
Scap = 25

apport_solaire=Scap * FSm * meteo[:,2]

besoinBrut = besoin_bat(Tconsigne,meteo[:,1],Rm,Ri,Rf) * agenda

besoin = besoinBrut - apport_solaire * agenda

plt.subplot(311)
plt.plot(Tsor_dro,label="Tsor_dro",color="red")
plt.plot(Tinj_dro,label="Tinj_dro",color="purple")
plt.legend()

plt.subplot(312)
plt.plot(Tinj_sto,label="Tinj_sto",color="orange")
plt.plot(Tsor_sto,label="Tsor_sto",color="blue")
plt.plot(meteo[i_summerStart:i_summerEnd,0],Tsable,label="Tsable",color="red")
plt.legend()

plt.subplot(313)
plt.plot(besoinBrut,label="besoin brut W")
plt.plot(besoin,label="besoin net W = besoin brut - apport solaire")
#plt.plot(meteo[:,2],label="apport solaires en W/m2")
plt.plot(apport_solaire,label="apport solaire en W")
plt.legend()

plt.show()


"""
# Pompe à Chaleur

COP=3
for i in range (0,meteo.shape[0]-1):
    
    Te2_sto[i]=Ts2_sto[i]-COP*Besoin[i]/((COP-1)*m2*Cp2)
    Ts1_sto[i]=Tsorties_stockage(m1,Cpf,Ts2_ech[i],Tsable[i],k_global)
    Ts2_sto[i]=Tsorties_stockage(m2,Cp2,Te2_sto[i],Tsable[i],k_global)
"""
    









