import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from dromosense.tools import *
from dromosense.constantes import rho_eau,Cpf,kelvin
from dromosense.sencity import *
import math
from datetime import datetime
from dateutil import tz
CET=tz.gettz('Europe/Paris')
import sys
from matplotlib.sankey import Sankey
"""
Fonctionnement en mode API
Le script peut être piloté en ligne de commande
Pour lancer le cas d'usage 0 avec ECS, il faut taper la commande suivante :
```
python couplage.py True 0
```
premier paramètre : True ou False pour activer ou non la production d'ECS
second paramètre : numéro du cas d'usage
Pour en plus enregistrer le graphe sortant dans le répertoire courant, 
il faut rajouter un troisième paramètre indiquant le nom à utiliser pour le graphique, sans préciser l'extension:
```
python couplage.py True 0 "SummerProdECS"
```
"""



def graphe(s,e):
    """
    courbes et graphiques représentant le système modélisé dans son ensemble
    """

    figure = plt.figure(figsize = (10, 10))
    matplotlib.rc('font', size=8)

    xrange = meteo[s:e,0]-s
    clearblue = "#8dd8ff"
    
    ax1 = plt.subplot(511)
    l1 = "couplage dromotherme/échangeur de séparation de réseau/stockage/PAC"
    legende = "RTD : rayonnement total sur dromotherme"
    legende = "{} / RDO : rayonnement sur dromotherme en opération".format(legende)
    legende = "{} / ECD : énergie captée par dromotherme".format(legende)
    if ECSupply:
        l2 = "température de consigne dans le bâtiment : {} °C - température ECS : {} °C".format(Tconsigne,Tballon)
    else:
        l2 = "température de consigne dans le bâtiment : {} °C".format(Tconsigne)
    if label:
        plt.title("{}\n{} -> {} - {} \n{}\n{}".format(l1,_s,_e,label,l2,legende))
    else:
        plt.title("{}\n{} -> {}\n{}\n{}".format(l1,_s,_e,l2,legende))
    
    ## graphe 1 - la route
    plt.ylabel('dromotherme\n{} m2'.format(int(Surface_dro)))
    
    ax1.plot(xrange, RSB.agenda_dro[s:e], color=clearblue, label="dro ON/OFF")
    ax1.legend(loc="upper left")
    
    ax11 = ax1.twinx()
    ax11.plot(xrange, RSB.Tsor_dro[s:e], label="Tsor_dro", color="red")
    ax11.plot(xrange, RSB.Tinj_dro[s:e], label="Tinj_dro", color="purple")
    ax11.legend(loc="upper right")  
    ax11.set_ylabel('°C')
    
    
    ## graphe 2 - l'énergie
    ax2 = plt.subplot(512, sharex=ax1)
    
    ax2.plot(xrange, Surface_dro*meteo[s:e,2], label="RTD", color="yellow")
    ax2.plot(xrange, Surface_dro*(RSB.agenda_dro*meteo[:,2])[s:e], label="RDO = {} KWh".format(int(Esolaire)), color="#f8cccc"); 
    ax2.plot(xrange, Pdro[s:e], label="ECD = {} kWh / rendement {} %".format(int(Edro), int(Taux)), color="green"); 
    plt.ylabel('Puissances \n W')
    ax2.legend()
    
    ## graphe 3 - le stockage
    ax3 = plt.subplot(513, sharex=ax1)
    plt.ylabel('stockage\n {} m^3'.format(int(V_sable)))
    
   # ax3.plot(xrange, RSB.diff[s:e], color=clearblue, label="derivée de Tsable en °C/s ou K/s")
    #ax3.legend(loc="upper left")
    
    ax31 = ax3.twinx()
    ax31.plot(xrange, RSB.Tinj_sto[s:e], label="Tinj_sto", color="orange")
    ax31.plot(xrange, RSB.Tsor_sto[s:e], label="Tsor_sto", color="blue")
    ax31.plot(xrange, RSB.Tsable[s:e], label="Tsable", color="red")
   
    
    ax31.legend(loc="upper right")
    ax31.set_ylabel('°C')
    
    ## graphe 4 - la PAC
    ax4 = plt.subplot(514, sharex=ax1)
    plt.ylabel('PAC')
    
    ax4.plot(xrange, RSB.agenda_pac[s:e], color="black", label="PAC ON/OFF")
    ax4.legend(loc="upper left")
    
    ax41=ax4.twinx()
    ax41.plot(xrange, RSB.Tinj_pac[s:e], label="Tinj_pac", color="red")
    ax41.plot(xrange, RSB.Tsor_pac[s:e], label="Tsor_pac", color="#7cb0ff")
    ax41.legend(loc="upper right")
    ax41.set_ylabel('°C')
    
    
    ## graphe 5 - le besoin du bâtiment
    ax5 = plt.subplot(515, sharex=ax1)
    plt.ylabel('Bâtiment \n W')
    plt.xlabel("Temps - 1 unité = {} s".format(step))
    
   
    
   
    ax5.plot(xrange, besoin_total[s:e], label="besoin {}  = {} kWh / Elec {} kWh".format(heating,int(conso_bat),int(Eelec)), color="orange")
    ax5.legend(loc="lower right")
    
    plt.show()
    
    if len(sys.argv) > 3:
        if sys.argv[3]:
            figure.savefig("{}.png".format(sys.argv[3]))
    
    

def ECSPower(min, max, size):
    """
    on modélise l'eau du réseau comme une fonction sinusoidale de période annuelle
    
    cette fonction est complètement calée sur un fichier météo qui commence au 1er janvier mais qui peut être pluriannuel
    
    min : température minimale d'injection de l'eau du réseau dans le ballon
    max : température maximale d'injection de l'eau du réseau dans le ballon
    """
    T_water=np.zeros(size)
    ## période
    w = 2*math.pi/npy
    for i in range(size):
        # numéro de step dans l'année
        siy = i - npy*(i//npy)     
        T_water[i]= 0.5 * ( (min-max)* math.cos(w*siy) + max + min )
    # le besoin s'entend pour une journée, ie 24*3600 secondes
    # il faut donc diviser par 24*3600 pour convertir de J à W, Cpf étant exprimée en J/kg/K
    return Volume_ballon*Npers*(Tballon-T_water)*Cpf/(24*3600)


"""
IMPORTATION DES DONNES METEOS
0 : temps exprime en heures
1 : temperature d'air (en deg Celsius)
2 : rayonnement global (en W/m2)
3 : rayonnement atmospherique (en W/m2)
4 : vitesse du vent (en m/s)
NOTA : principe d'un échantillonnage temporel à l'année = il faut donner à ce script un fichier meteo annuel
 
la variable npy représente le nombre de points dans une année complète
dans le cas présent, l'intervalle de temps est l'heure. 
Si on choisit un autre pas, il faut changer la variable step, exprimée en secondes, manuellement !!!!
Dansc ce cas, si on désigne par dt1 le nouveau pas de temps, il va falloir multiplier le pas d'espace dx par dt1/3600 
pour respecter la condition de stabilité du schema numérique utilisé au niveau du domotherm.
"""

#meteo = np.loadtxt('../../datas/meteo_Bourget_du_lac.txt') 
meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')
npy = meteo.shape[0]
##meteo=np.concatenate((meteo,meteo,meteo,meteo,meteo,meteo,meteo,meteo,meteo,meteo))
year=3
meteo=np.concatenate((meteo,meteo,meteo))# La simulation sur plus de deux années permet de faire converger les solutions
for i in range(year):
    meteo[i*npy:(i+1)*npy,0]=meteo[i*npy:(i+1)*npy,0]+i*npy

print(meteo.shape)
f2 = 1000.0*1.1*(0.0036*meteo[:,4]+0.00423)
f1 = (1.0-albedo)*meteo[:,2] + meteo[:,3] + f2*(meteo[:,1]+kelvin)

# longueur et largeur de l'échangeur, exprimées en m
lincha=7.5
larcha=4
Surface_dro=larcha*lincha
dx=0.75
# pas de temps en secondes pour la discrétisation temporelle
step=3600

"""
débit dans le dromotherme
qdro_u est le débit unitaire dans le dromotherme en m3/s :
- pour un mètre linéaire de chaussée selon le profil en long 
- pour une largeur de chaussée de 4 mètres selon le profil en travers 
"""
qdro_u = 0.035/step # m3/s
qdro=qdro_u*lincha
start = 1483232400
summerStart=1493600400 # 1er mai 
summerEnd = 1506819600 # 30 septembre


# instanciation d'un dromotherme 1D - input.txt contient les paramètres calant ce modèle sur le modèle 2D
dromo=OneDModel('input.txt',step,meteo.shape[0],larcha,dx)
dromo.f1 = f1
dromo.f2 = f2
#dromo.T[0,:,:] = np.ones((dromo.T.shape[1],dromo.T.shape[2]))*10+kelvin
# très provisoire, il faudrait discuter de celà
dromo.T = np.ones((dromo.T.shape[0],dromo.T.shape[1],dromo.T.shape[2]))*10+kelvin

# route stockage bâtiment
RSB=SenCityOne(meteo.shape[0],step)

"""
massif de stockage
"""
m_sable=70200.0 # masse de sable humide (sable sec + eau ) en kg
Cp_sable=1470.0 # capacité calorifique massique du sable humide en J/Kg.K
V_sable=45 # volume du stockage (en m^3)
"""
paramètres définissant le système géothermique équipant le stockage
"""
R1=0.0102 # Rayon intérieur du tube
R2=0.0125 # Rayon extérieur
L_tube=5 # Longueur tube en m
N_tube=16  # Nombre de tubes
# conductivités en W/(K.m)
lambda_pac=15.8
lambda_tube=1.32
Nu=4.36 # Nombre de Nusselt de l'écoulement du fluide à l'intérieur des tubes
Rcond=np.log(R2/R1)/(2*math.pi*N_tube*L_tube*lambda_tube)
Rconv=1/(math.pi*Nu*N_tube*L_tube*lambda_pac)
# k exprimé en W/K
k=1/(Rcond+Rconv)

# efficacité de l'échangeur
eff = 0.8
# débit dans la partie de l'échangeur côté "stockage"
qsto = 1.2*qdro

# rho_eau en provenance du fichier constantes est exprimée en kg/m3
# Cpf en provenance du fichier constantes est exprimée en J/(kg.K)
mdro = qdro * rho_eau
msto = qsto * rho_eau
cpdro = Cpf
cpsto = Cpf

coeff = (mdro * cpdro) / (msto * cpsto)

"""
besoin en chauffage
"""
Tconsigne=19
# toutes les résistances sont exprimées en K/W
Rm=8.24E-02 # Résistance thermique des murs
Ri=1.43E-03 # Résistance superficielle intérieure
Rf=3.4E-02 # Résistance due aux infiltrations+vitre et au renouvellement d'air
# facteur solaire pour des matériaux plein de type murs/toits
FSm=0.0048
"""
hypothèse : la maison fait 4 mètres de large sur 5 mètres de long
vu du ciel c'est donc un rectangle de 20m2 de surface
"""
Scap = 20

"""
PAC
"""
COP=3
rho_pac=1040 # en kg/m^3
mpac=rho_pac*0.035*4/3600
cpac=3942.0  #' en J/kg.K

Tmoy=11
Tamp=9.1
lambda_sable=1.7
rho_sable=1700
e_iso=0.2
SL_iso=40.5
SB_iso=20
lambda_iso=0.038


RSB.set(eff,k,coeff,msto,cpsto,m_sable,Cp_sable,mpac,cpac)
RSB.setPertes(Tmoy,Tamp,lambda_sable,rho_sable,e_iso,SL_iso,SB_iso,lambda_iso)


"""
on calcule les index permettant de boucler sur l'été, ainsi que le besoin net du bâtiment et la puissance géothermique à développer
"""
i_summerStart=(summerStart-start)//step

i_summerEnd=i_summerStart+(summerEnd-summerStart)//step

apport_solaire = Scap * FSm * meteo[:,2]

nbpts=meteo.shape[0]

#schedule_bureau=np.array([[8,17],[8,17],[8,17],[8,17],[8,17],[-1,-1],[-1,-1]])

#agenda_bureau=basicAgenda(nbpts,step,start,summerStart,summerEnd,schedule=schedule_bureau)

#schedule_domestique=np.array([[18,23],[18,23],[18,23],[18,23],[18,23],[7,23],[7,23]])

#agenda_domestique=basicAgenda(nbpts,step,start,summerStart,summerEnd,schedule=schedule_domestique)

besoinBrut = besoin_bat(Tconsigne,meteo[:,1],Rm,Ri,Rf)

besoin_chauffage = besoinBrut - apport_solaire

"""
on efface les besoins de chauffage sur les étés
"""
for i in range(year):
#besoin_chauffage[i_summerStart:i_summerEnd] = np.zeros(i_summerEnd-i_summerStart)

    besoin_chauffage[i_summerStart+i*npy:i_summerEnd+i*npy]=np.zeros(i_summerEnd-i_summerStart)

"""
Besoin en ECS
On considère que le besoin en ECS d'une personne seule oscille entre 30 et 40L.
Nous prendrons donc 35L comme besoin journalier.
La température de l'eau chaude stockée dans le ballon est de 60°C.
L'eau entre dans le ballon à une température de 10°C en huvers et 16°C en été
"""
Tballon=60
Tentree_ete=16  # Tmax
Tentree_hiver= 10  # Tmin
Volume_ballon=35 # 35L/jour en moyenne pour une personne
Npers=6

ECS = ECSPower(Tentree_hiver, Tentree_ete, meteo.shape[0])

"""
*************************************
*************************************
!!!!!!!!!!! SWITCHES !!!!!!!!!!!
ON PEUT ACTUELLEMENT SIMULER 8 CAS DE FIGURES
- choix du usecase
- fait-on de l'ECS ou pas ?
exemples :
1) usecase=0 avec ECSupply=False => recharge du stock été sans consommation aucune
2) usecase=0 avec ECSupply=True => consommation estivale d'ECS
3) usecase=1 avec ECSSupply=False => dromotherme été pour recharge du stock puis utilisation pour chauffage sur début hiver
4) usecase=1 avec ECSSupply=True => dromotherme été pour recharge du stock puis utilisation pour chauffage+ECS sur début hiver
5) usecase=2 avec ECSSupply=False => dromotherme toute l'année et utilisation pour chauffage sur hiver entier
6) usecase=2 avec ECSSupply=True => dromotherme toute l'année et utilisation pour chauffage+ECS sur hiver entier
7) usecase=3 avec ECSSupply=False => dromotherme été + hiver si rayonnement au dessus d'un seuil et utilisation pour chauffage sur hiver entier
8) usecase=3 avec ECSSupply=True => dromotherme été + hiver si rayonnement au dessus d'un seuil et utilisation pour chauffage+ECS sur hiver entier
"""
if len(sys.argv) > 2:
    ECSupply= (sys.argv[1] == "True")
    usecase= int(sys.argv[2])
else:    
    ECSupply=True
    usecase=0
"""
*************************************
*************************************
"""
            
if ECSupply:
    heating = "chauffage+ECS"
else:
    heating = "chauffage"
            
besoin_total = (besoin_chauffage + ECSupply * ECS)

RSB.Pgeo = (COP-1) * besoin_total / COP

label = ""


"""
Management des usecases
1) fixer les indices simStart et simEnd pour définir la fenêtre de simulation
2) définir les agendas agenda_dro et agenda_pac
3) réajuster éventuellement la chaine heating si l'on souhaite autre chose que les valeurs par défaut définies plus haut
4) réajuster éventuellement la chaine label qui est vide si on souhaite afficher une singularité dans le graphique
POUR CREER UN NOUVEAU USECASE, IL FAUT DONC A MINIMA DEFINIR simStart, simEnd, agenda_dro et agenda_pac
"""

if usecase == 0:
    # ECS only pendant l'été
    simStart = i_summerStart
    simEnd = i_summerEnd+365*24*(year-1)
    k1=simStart+365*24*(year-2)
    k2=simEnd-365*24*(year-2)
    RSB.agenda_dro[simStart:simEnd]=np.ones(simEnd-simStart)
    RSB.agenda_pac[simStart:simEnd]=np.ones(simEnd-simStart)
    heating = "ECS"

if usecase == 1:
    # dromotherme durant l'été
    simStart = i_summerStart
    simEnd=i_summerStart+365*24*(year-1)
    k1=simStart+365*24*(year-2)
    k2=simEnd
    RSB.agenda_dro[simStart:i_summerEnd]=np.ones(i_summerEnd-simStart)
    RSB.agenda_dro[simStart+365*24*1:i_summerEnd+365*24*1]=np.ones(i_summerEnd-simStart)
    RSB.agenda_pac[simStart:simEnd]=np.ones(simEnd-simStart)

if usecase == 2:
    # simulation annuelle
    simStart = i_summerStart
    simEnd=i_summerStart+365*24*(year-1)
    k1=simStart+365*24*(year-2)
    k2=simEnd
    RSB.agenda_dro[simStart:simEnd]=np.ones(simEnd-simStart)
    RSB.agenda_pac[simStart:simEnd]=np.ones(simEnd-simStart)

if usecase == 3:
    # simulation annuelle
    # dromotherme l'été et par intermittence l'hiver quant le rayonnement global est au dessus de 250 W/m2
    simStart = i_summerStart
    simEnd=i_summerStart+365*24*(year-1)
    k1=simStart+365*24*(year-2)
    k2=simEnd
    RSB.agenda_dro[simStart:i_summerEnd]=np.ones(i_summerEnd-simStart)
    RSB.agenda_pac[simStart:simEnd]=np.ones(simEnd-simStart)
    level=250
    label = "hiver si ray.>250 W/m2 : dromo on".format(level) 
    for i in range(i_summerEnd,simEnd):
        if meteo[i,2] >= level:
            RSB.agenda_dro[i]=1
 
    
if usecase == 4:
    simStart = i_summerStart
    simEnd=i_summerStart+365*24*(year-1)
    k1=simStart+365*24*(year-2)
    k2=simEnd
    RSB.agenda_dro[simStart:simEnd]=np.ones(simEnd-simStart)
    RSB.agenda_pac[simStart:simEnd]=np.zeros(simEnd-simStart)
    heating = "Pas de besoin"    
    
    
for i in range(simStart,simEnd):
    if RSB.Pgeo[i]==0:
        RSB.agenda_pac[i]=0            


"""
timestamps de la fenêtre de simulation + expressions compréhensibles pour humains
"""
tss= start + simStart * step
tse= start + simEnd * step
_s=tsToHuman(tss)
_e=tsToHuman(tse)


for i in range (int(simStart),int(simEnd)):
   RSB.SystemLoop(i,qdro_u,dromo)

"""
BILAN ENERGETIQUE
dromotherme
"""

Pdro = mdro * cpdro * RSB.agenda_dro * (RSB.Tsor_dro - RSB.Tinj_dro )
# toutes valeurs en kWh/m^2
Edro=np.sum(Pdro[k1:k2])*step/(3600*1000)
Esolaire=(1-albedo)*np.sum(step*RSB.agenda_dro[k1:k2]*meteo[k1:k2,2]*Surface_dro)/(3600*1000)
# taux de récupération en %
Taux=Edro*100/Esolaire

"""
bâtiment
"""
conso_bat=np.sum(RSB.agenda_pac[k1:k2]*besoin_total[k1:k2])*step/(3600*1000)
Eelec=conso_bat/COP
Eprimaire=2.58*Eelec

"""
LES PERTES
Version provisoire du calcul des pertes. Cette partie sera "transférée" vers un des modules
"""

Pertes_rad=np.zeros(meteo.shape[0])# Les pertes radiatives instantanées
Pertes_conv=np.zeros(meteo.shape[0])# Les pertes convectives instantanées

for i in range (k1,k2):
    Pertes_rad[i]=RSB.agenda_dro[i]*(epsilon*sigma*(dromo.T[i,0,-1])**4-meteo[i,3])
    Pertes_conv[i]=RSB.agenda_dro[i]*(1000.0*1.1*(0.0036*meteo[i,4]+0.00423))*(dromo.T[i,0,-1]-kelvin-meteo[i,1])

Prad=Surface_dro*np.sum(Pertes_rad)*step/(3600*1000) # Pertes radiatives totales
Pconv=Surface_dro*np.sum(Pertes_conv)*step/(1000*3600) # Pertes convectives totales
Psto=np.sum(RSB.pertes[k1:k2])*step/(1000*3600) # Pertes totales au niveau du stockage
Egeo=np.sum(RSB.agenda_pac[k1:k2]*RSB.Pgeo[k1:k2])*step/(1000*3600) # Energie géothermique 
Echauff=np.sum(RSB.agenda_pac[k1:k2]*besoin_chauffage[k1:k2])*step/(1000*3600) # Consommation du chuaffage
print("Esolaire=",Esolaire)
print("Prad=", Prad)
print("Edro=", Edro)
print("Pconv=", Pconv)
print("Psto", Psto)

integralite=False # booléen donnant l'intégralité ou non du graphique
if integralite:
    k1=simStart
    k2=simEnd

graphe(k1,k2)    
    