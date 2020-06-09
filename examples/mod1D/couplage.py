import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from dromosense.tools import *
from dromosense.constantes import *
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

#meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')
meteo = np.loadtxt('meteo2.txt')
print(meteo.shape)
#T = np.loadtxt('T1d.txt')#+kelvin
#print(T.shape)
input("press any key")


def F(y,t):
    """
    y : Tsable = Tstockage

    t : time index

    result : dy/dt = dTsable/dt = dTstockage/dt
    """
    #i = int(t/3600)
    i = int(t)
    if verbose:
        print("we have t={} and y={}".format(i,y))

    #Tsor_sto[i] = ( k * y / (msto * cpsto + k/2) - coeff * eff * Tsor_dro[i] ) / a

    #Tsor_sto[i] = ( k * y - B * Tsor_dro[i] ) / A
        
    dromo.iterate(i,Tinj_dro[i])
    Tsor_dro[i]=dromo.T[i,1,-1]-kelvin    
    Tsor_sto[i] = ( k * y + B * Tsor_dro[i] ) / A

    Tinj_sto[i] = Tsor_sto[i] + coeff * eff * (Tsor_dro[i] - Tsor_sto[i])

    Tinj_dro[i] = Tsor_dro[i] - eff * (Tsor_dro[i] - Tsor_sto[i])

    Tinj_pac[i]=y-C*Pgeo[i]/k

    Tsor_pac[i]=Tinj_pac[i]-Pgeo[i]/(mpac*cpac)
     
    #der=msto * cpsto * (Tinj_sto[i] - Tsor_sto[i]) / (m_sable * Cp_sable)
    # agenda n'est pas utile là.....Pgeo intègre déjà les effets de l'agenda vu sa construction.....
    #der=(msto * cpsto * (Tinj_sto[i] - Tsor_sto[i]) + mpac*cpac*(Tsor_pac[i]-Tinj_pac[i])) / (m_sable * Cp_sable)
    
    # tu n'as pas besoin de refaire ce calcul de mpac*cpac*(Tsor_pac[i]-Tinj_pac[i])
    # normalement vu la ligne 51, mpac*cpac*(Tsor_pac[i]-Tinj_pac[i]) vaut exactement Pgeo[i]
    # pourquoi n'obtient-on pas les mêmes résultats quant on utilise la ligne 59 en lieu et place de la ligne 55 ?
    der = ( msto * cpsto * (Tinj_sto[i] - Tsor_sto[i]) + Pgeo[i] ) / (m_sable * Cp_sable)

    if verbose:
        print("dTsable/dt is {}".format(der))

    return der


# température d'entrée et de sortie du fluide dans le stockage
Tinj_sto=np.zeros(meteo.shape[0]+1)
# température de sortie du fluide après transit dans le stockage
Tsor_sto=np.zeros(meteo.shape[0])
# température d'entrée du fluide géothermique dans le stockage (sortie de la PAC)
Tsor_pac=np.zeros(meteo.shape[0])
# température de sortie du fluide géothermique dans le stockage (entrée  de la PAC)
Tinj_pac=np.zeros(meteo.shape[0])
# température d'injection du fluide dans le dromotherme
Tinj_dro=10*np.ones(meteo.shape[0])
# température de sortie de fluide après collecte de chaleur dans le dromotherme
#Tsor_dro=T[:,1]
Tsor_dro=np.zeros(meteo.shape[0])


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

#B = (msto * cpsto +k/2) * coeff * eff
B = (msto * cpsto -k/2) * coeff * eff
#A = k - B
A = k + B

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
#Tsable = odeint(F,10,meteo[i_summerStart:i_summerEnd,0])

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
Ci=18.407 # Capacité thermique de l'air (Wh/K)
Cm=2636 # Capacité thermique des murs
Sm=49.5
FSm=0.0048
Sv=3.5
FSv=0.8
"""
hypothèse : la maison fait 5 mètres de large sur 5 mètres de long
vu du ciel c'est donc un carré de 25m2 de surface
"""
Scap = 20

#apport_solaire=Scap * FSm * meteo[:,2]
apport_solaire=Scap * FSm * meteo[:,5]
besoinBrut = besoin_bat(Tconsigne,meteo[:,1],Rm,Ri,Rf) * agenda

besoin = besoinBrut - apport_solaire * agenda


# PAC
COP=3
#Pgeo=COP*besoin/(COP-1)
Pgeo=(COP-1)*besoin/COP
mpac=msto

cpac=4180.0
C=1-k/(2*mpac*cpac)


"""
IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
0 : temps exprime en secondes depuis le 01/06 00:00
1 : temperature d'air (en deg Celsius)
2 : temperature du point de rosee (pas utile)
3 : nature des precipitations (pas utile)
4 : vitesse du vent (en m/s)
5 : rayonnement global (en W/m2)
6 : rayonnement atmospherique (en W/m2)
"""

f2 = 1000.0*1.1*(0.0036*meteo[:,4]+0.00423)
f1 = (1.0-albedo)*meteo[:,5] + meteo[:,6] + f2*(meteo[:,1]+kelvin)
t = meteo[:,0]
nt = t.size # le nombre de points dans la discretisation temporelle
dt = 3600.0

"""
l'initialisation complète de l'objet se fait de la manière suivante :
dromo=OneDModel('input.txt',dt,nt,L,dx,qf)
avec :
- dt : pas de temps en secondes
- nt : nombre de points dans la discrétisation temporelle
- L : largeur de la chaussée en m
- dx : pas d'espace en m
- qf : débit en m3
le fluide est de l'eau
la classe OneDModel va chercher les valeurs de Cpf et rho_eau dans le fichier des constantes
içi on ne donne pas L, dx ni qf car la classe dispose de valeurs par défaut (4,0.75,0/035/3600)
mais dans la pratique, il faut fournir les bonnes valeurs

ensuite :
- il faut affecter les vecteurs dromo.f1 et dromo.f2 issus de l'exploitation du monitoring :-)
- il faut fixer la condition initiale pour dromo.T[0,:,:] - cf plus loin en utilisant T2d[0,1:]

Nota : les variables de classe sont publiques en python, contrairement au C++
on y accède via nom_instance.nom_variable
içi l'instance est dromo
"""
# instanciation.....
dromo=OneDModel('input.txt',dt,nt)
#dromo.showVectors()
dromo.f1=f1
dromo.f2=f2
input("press any key")


"""
IMPORTATION DES TEMPERATURES DU MODELE 2D
dans l'ordre
colonne 0 : temps
colonne 1 : température couche de surface
colonne 2 : température couche drainante
colonne 3 : température couche de base
colonne 4 : température couche fictive ?
colonne 5 : température massif
"""

T2d = kelvin + np.loadtxt('T2d2.txt')

L = 4.0 # Largeur de chaussee en m
dx = 0.75 # pas d'espace en m
qf = 0.035/3600.0         # debit volumique du fluide (m^3/s)
Cf = Cpf*rho_eau #4200000.0 # capacite calorifique volumique de l'eau (J/(m3.K))
phi = 0.0     #( porosite de la couche drainante)

nx = int(L/dx) # Nombre de blocs

# Conditions initiales
for i in range(nx) :
   # T[0,:,i] = T2d[0,1:]
    dromo.T[0,:,i] = T2d[0,1:]

# Conditions aux limites
#Tinj = 10.0 + kelvin

#for n in range(1,nt):
    #dromo.iterate(n,Tinj)




#Température du stockage/sable
#pourquoi le stock serait-il à 10°C au milieu de l'hiver??? c'est au sortir de l'hiver quon pense qu'il sera à 10
Tsable = odeint(F,10,meteo[:,0])
Pgeo2=mpac*cpac*(Tinj_pac-Tsor_pac)
Tsor_pac_wastewater=10-Pgeo/(mpac*cpac)




figure = plt.figure(figsize = (10, 10))
matplotlib.rc('font', size=8)

ax1 = plt.subplot(411)
ax1.plot(Tsor_dro,label="Tsor_dro",color="red")
ax1.plot(Tinj_dro,label="Tinj_dro",color="purple")
ax1.legend()

ax2 = plt.subplot(412, sharex=ax1)
ax2.plot(Tinj_sto,label="Tinj_sto",color="orange")
ax2.plot(meteo[:,0],Tsable,label="Tsable",color="k")
ax2.plot(Tsor_sto,label="Tsor_sto",color="blue")
ax2.legend()

ax3 = plt.subplot(413, sharex=ax1)
#ax3.plot(Tsor_pac_wastewater,label="Tsor_pac_wastewater10°C",color="blue")
ax3.plot(Tinj_pac,label="Tinj_pac",color="red")
ax3.plot(Tsor_pac,label="Tsor_pac",color="#7cb0ff")
ax3.plot(Pgeo-Pgeo2,label="Ecart",color="k")
#plt.plot(meteo[i_summerStart:i_summerEnd,0],Tsable,label="Tsable",color="red")
ax3.legend()

ax4 = plt.subplot(414, sharex=ax1)
ax4.plot(besoinBrut,label="bes.brut W",color="red")
ax4.plot(besoin,label="bes.net W = bes.brut - app.sol",color="orange")
#ax4.plot(besoin_surfacique,label="besoin net par unité de surface W/m^2",color="orange")
#ax4.plot(meteo[:,2],label="app.sol W/m2")
ax4.plot(apport_solaire,label="app.sol en W",color="yellow")
#ax4.plot(Pdro,label="Densité de lux du dromotherm en W/m^2",color="green")
ax4.legend()

plt.show()


plt.subplot(111)
plt.title("modèle 1D vs classe")
plt.ylabel("°C")
plt.xlabel("heures")
#plt.plot(t/3600,T[:,1,-1]-kelvin,label="1D model")
plt.plot(t/3600,dromo.T[:,1,-1]-kelvin,label="classe")
plt.legend(loc="upper right")
plt.show()








