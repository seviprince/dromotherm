import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from dromosense.tools import *
from dromosense.constantes import rho_eau,Cpf,kelvin


from dromosense.sencity2 import *
from dromosense.fstorage import *
import math


from datetime import datetime
from dateutil import tz
CET=tz.gettz('Europe/Paris')
import sys
from matplotlib.sankey import Sankey

import datetime as dt
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

    ## graphe 1 - la route
    plt.ylabel('dromotherm')
    
    ax1.plot(xrange, RSB.agenda_dro[s:e], color=clearblue, label="dro ON/OFF")
    ax1.legend(loc="lower left")
    
    ax11 = ax1.twinx()
   # ax11.plot(xrange, Tdro_expe, label="T drainant expé", color="black")
    ax11.plot(xrange, RSB.Tsor_dro[s:e], label="Tsor_dro", color="red")
    ax11.plot(xrange, RSB.Tinj_dro[s:e], label="Tinj_dro", color="purple")
    ax11.legend(loc="upper right")  
    ax11.set_ylabel('°C')
    
    
 
    
    ## graphe 2 - Le stocakge/ Températures couches non englacées
    ax2 = plt.subplot(512, sharex=ax1)
    labels=["Tl1","Tl2","Tl3","Tl4","Tl5"]
    colors=["r","k","r--","k--","b"]
    plt.ylabel('couche stockage \n °C')
    plt.ylabel("Température (°C)")
    for j in range(0,nc):
        ax2.plot(RSB.Tl[s:e,j],colors[j],label=labels[j])

    plt.xlabel("Time - 1 unit = {} s".format(step)) 
    ax2.legend(loc="lower right")
    
    ## graphe 3 - Le stocakge/ Températures aux interfaces des couches
    ax3 = plt.subplot(513, sharex=ax1)
    labels=["T1/2","T2/3","T3/4","T4/5"]
    plt.ylabel("Température (°C)")
    for j in range(2,6):
        ax3.plot(RSB.T_inter[s:e,j],colors[j-2],label=labels[j-2])
    plt.xlabel("Time - 1 unit = {} s".format(step)) 
    ax3.legend(loc="lower right")

    ## graphe 4 - le stockage/ Température fluide_dromotherm
    ax4 = plt.subplot(514, sharex=ax1)
    plt.ylabel('stockage \n °C')
   # ax2.plot(xrange, RSB.Tmoy_inj_sto[s:e], label="Tinj storage", color="red")   
   # ax2.plot(xrange, RSB.Tmoy_sor_sto[s:e], label="Tsor storage", color="blue") 
   # ax2.plot(xrange, np.mean(RSB.Tl[s:e,:]), label="T storage", color="k") 
    ax4.plot(xrange, RSB.Tmoy_inj_sto[s:e], label="Tinj sto", color="red")   
    ax4.plot(xrange, RSB.Tmoy_sor_sto[s:e], label="Tsor sto", color="blue")          
    ax4.legend(loc="upper right")

    ## graphe 5 -la PAC
    ax5 = plt.subplot(515, sharex=ax1)
    plt.ylabel('stockage \n °C')
   # ax2.plot(xrange, RSB.Tmoy_inj_sto[s:e], label="Tinj storage", color="red")   
   # ax2.plot(xrange, RSB.Tmoy_sor_sto[s:e], label="Tsor storage", color="blue") 
   # ax2.plot(xrange, np.mean(RSB.Tl[s:e,:]), label="T storage", color="k") 
    ax5.plot(xrange, RSB.Tmoy_inj_pac[s:e], label="Tinj pac", color="red")   
    ax5.plot(xrange, RSB.Tmoy_sor_pac[s:e], label="Tsor pac", color="blue")          
    ax5.legend(loc="upper right")
    
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
from PyFina import getMeta, PyFina


# Feeds ID
feeds = {
    "temp_ext": 13,
    "netRad":14,
    "ray_atmo":81,
    "windSpeed": 84,
    }
step=3600
# Path of the PyFina File
dir = "D:/Utilisateurs/sevif/Documents/GitHub/Tuto_Pierre/data/emoncms-backup-2023-03-09/phpfina/"



##meteo=np.concatenate((meteo,meteo,meteo,meteo,meteo,meteo,meteo,meteo,meteo,meteo))
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
nc=5
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
#RSB=SenCityOne(meteo.shape[0],step)

RSB=SenCityTwo(step,meteo.shape[0],nc)
RSB.Text=meteo[:,1]
y1=RSB.Text
"""
massif de stockage
"""

V_sable=45 # volume du stockage (en m^3)
"""
paramètres définissant le système géothermique équipant le stockage
"""
R1=0.0102 # Rayon intérieur du tube
R2=0.0125 # Rayon extérieur
L_tube=100 # Longueur tube en m
N_tube=1  # Nombre de tubes
# conductivités en W/(K.m)
lambda_pac=0.4
lambda_tube=1.32
Nu=4.36 # Nombre de Nusselt de l'écoulement du fluide à l'intérieur des tubes
Rcond=np.log(R2/R1)/(2*math.pi*N_tube*L_tube*lambda_tube)
Rconv=1/(math.pi*Nu*N_tube*L_tube*lambda_pac)
Rthf=Rcond+Rconv
# k exprimé en W/K
k=1/(Rcond+Rconv)

# efficacité de l'échangeur
eff = 0.8
# débit dans la partie de l'échangeur côté "stockage"

qsto=0.3/3600
# rho_eau en provenance du fichier constantes est exprimée en kg/m3
# Cpf en provenance du fichier constantes est exprimée en J/(kg.K)
mdro = qdro * rho_eau
msto = qsto * rho_eau
cpdro = Cpf




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
#mpac=rho_pac*0.035*4/3600
mpac=rho_pac*0.75/3600
cpac=3942.0  #' en J/kg.K
cpsto = cpac
Tmoy=11
Tamp=9.1
lambda_sable=1.7
rho_sable=1700
e_iso=0.2
SL_iso=40.5
SB_iso=20
lambda_iso=0.038

#paramètres thermo-physiques des différents élements du stockage


"""
Le béton
"""
beton = {"e":0.2, "lambda":2.0, "rho":2500.0, "cp":1000.0}
ath_be=beton["lambda"]/(beton["rho"]*beton["cp"]) # diffusivité thermique du béton

"""
Isolant
"""
isolant={"e":0.2, "lambda":0.038, "rho":50.0, "cp":1450.0}
ath_iso=isolant["lambda"]/(isolant["rho"]*isolant["cp"]) # diffusivité thermique du de l'isolant


"""
Couche de sable + eau
"""
sable={"lambda":1.7,"rho":1700.0,"cp":1470.0}

eau={"lambda":0.56,"rho":1000.0,"cp":4181.0}

glace={"lambda":2.24,"rho":910.0,"cp":2100}
L=333000 #chaleur latente de fusion de l'eau à 0°c en J/kg


"""
paramètres thermo-physiques des couches

lambda,rho, Cp, m,et a signifient respectivement: conductivité thermique, masse volumique,capacité thermique massique
masse et diffusivité thermique.

Indice s : solide (sable+glace)

Indice l : liquide (sable+eau liquide)

"""
phi=0.27 # taux de vide dans chaque couche

lambda_s=lambda_eq(glace["lambda"],sable["lambda"],phi)
lambda_l=lambda_eq(eau["lambda"],sable["lambda"],phi)

rho_s=(1-phi)*sable["rho"]+phi*glace["rho"]
rho_l=(1-phi)*sable["rho"]+phi*eau["rho"]

Cp_s=(1-phi)*sable["cp"]+phi*glace["cp"]
Cp_l=(1-phi)*sable["cp"]+phi*eau["cp"]


ath_s=lambda_s/(rho_s*Cp_s)
ath_l=lambda_l/(rho_l*Cp_l)


#Données sur le stockage
nc=5 #nombre de couches


lambda_be=2.0
lambda_iso=0.038
cste=step/(phi*glace["rho"]*L)
coeff = (mdro * cpdro) / (msto * cpsto)
Sb=25
RSB.set(eff,k,coeff,cste,lambda_be,ath_be,lambda_iso,ath_iso,lambda_l,ath_l,rho_l,Cp_l,lambda_s,ath_s,rho_s,msto,cpsto,mpac,cpac,Rthf,Sb)



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
Npers=4

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
    ECSupply=False
    usecase=5
"""
*************************************
*************************************
"""
            
if ECSupply:
    heating = "chauffage+ECS"
else:
    heating = "chauffage"
            


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
    RSB.agenda_dro[simStart:simEnd]=0*np.ones(simEnd-simStart)
    RSB.agenda_pac[simStart:simEnd]=0*np.ones(simEnd-simStart)

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


if usecase == 5:

    simStart =1385+365*24# 27 février 
    simEnd=1608+365*24 # 09 mars
    #simStart =100000
    #simEnd=simStart+246
    k1=simStart
    k2=simEnd
    RSB.agenda_dro[simStart:simEnd]=0*np.ones(simEnd-simStart)
    RSB.agenda_pac[simStart:simEnd]=1*np.ones(simEnd-simStart)
    heating = "heating"        
    
for i in range(simStart,simEnd):
    if RSB.Pgeo[i]==0:
        RSB.agenda_pac[i]=0       

feeds = {
    "PaEnStPT": 100,
    "PaSoStPT": 102,
    
    "PaSoBaPT": 110,
    "PaEnBaPT": 107,
    
    "BaEnVcPT":109,
    "BaSoVcPT":108,
    
    "ChDrAmTh1":30,"ChDrAmTh2":32,
    "ChDrCeTh1":72, "ChDrCeTh2":67,
    "ChDrAvTh2":73, "ChDrAvTh4":74,
    
    "StC1EsOu1":48,"StC1EsOu2":46,"StC1EsOu3":47,"StC1EsOu4":49,"StC1EsOu5":50,
    "StC2EsOu1":51,"StC2EsOu2":54,"StC2EsOu3":52,"StC2EsOu4":56,"StC2EsOu5":55,    
    "StC3EsOu1":59,"StC3EsOu2":61,"StC3EsOu3":57,"StC3EsOu4":60,"StC3EsOu5":58,   
    "StC4EsOu1":45,"StC4EsOu2":42,"StC4EsOu3":52,"StC4EsOu4":43,"StC4EsOu5":44,      
    "StC5EsOu1":64,"StC5EsOu2":63,"StC5EsOu3":62,"StC5EsOu4":65,"StC5EsOu5":66,
    "P1":90
    }
# heures de début et fin des données 

localstart_string="2023-02-27 17:00:00" # heure de début 
localend_string="2023-03-09 23:00:00" # heure de fin
timestamp_start =int( dt.datetime.timestamp(dt.datetime.strptime(localstart_string, '%Y-%m-%d %H:%M:%S')))
timestamp_end =int( dt.datetime.timestamp(dt.datetime.strptime(localend_string, '%Y-%m-%d %H:%M:%S')))
nbpts=(timestamp_end -timestamp_start)//step
npt=k2-k1
#Tdro_expe=np.zeros(npy)
Tdro_expe=PyFina(feeds["ChDrAmTh1"],dir,timestamp_start ,step,npt)+PyFina(feeds["ChDrAmTh2"],dir,timestamp_start ,step,npt) \
    +PyFina(feeds["ChDrCeTh1"],dir,timestamp_start ,step,npt)+PyFina(feeds["ChDrCeTh2"],dir,timestamp_start ,step,npt) \
    +PyFina(feeds["ChDrAvTh2"],dir,timestamp_start ,step,npt)+PyFina(feeds["ChDrAvTh4"],dir,timestamp_start ,step,npt)    

Tdro_expe=Tdro_expe/6

Tsto=PyFina(feeds["StC1EsOu1"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC1EsOu2"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC1EsOu3"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC1EsOu4"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC1EsOu5"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC2EsOu1"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC2EsOu2"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC2EsOu3"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC2EsOu4"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC2EsOu5"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC3EsOu1"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC3EsOu2"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC3EsOu3"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC3EsOu4"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC3EsOu5"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC4EsOu1"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC4EsOu2"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC4EsOu3"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC4EsOu4"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC4EsOu5"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC5EsOu1"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC5EsOu2"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC5EsOu3"],dir,timestamp_start ,step,npt)+PyFina(feeds["StC5EsOu4"],dir,timestamp_start ,step,npt)\
    +PyFina(feeds["StC5EsOu5"],dir,timestamp_start ,step,npt)\
  
Tsto=Tsto/25  
Tpac_sor_sto=PyFina(feeds["PaSoStPT"],dir,timestamp_start ,step,npt)
Tpac_ent_sto=PyFina(feeds["PaEnStPT"],dir,timestamp_start ,step,npt)
   
Tpac_sor_ba=PyFina(feeds["PaEnBaPT"],dir,timestamp_start ,step,npt)
Tpac_ent_ba=PyFina(feeds["PaSoBaPT"],dir,timestamp_start ,step,npt)

Tba_sor_vc=PyFina(feeds["BaSoVcPT"],dir,timestamp_start ,step,npt)
Tba_ent_vc=PyFina(feeds["BaEnVcPT"],dir,timestamp_start ,step,npt)

Pelec=PyFina(feeds["P1"],dir,timestamp_start ,step,npt)/1000
#----------------------------------------------------------
rho=1000 # kg/m3
msto=rho*0.67/3600 #(0.67m3/heure)
mchaud=rho*0.55/3600 # 0.55m3/heure
mvc=rho_pac*0.4/3600 # 0.4m3/heure
Cpg=3942.0
Cpe=4180
Pelec2=np.zeros(npy)    
#----------------------------------------------------------


Pchaud=np.zeros(npt)
Pvc=np.zeros(npt)
COP=np.zeros(npt)
COP2=np.zeros(npt) 
       

  
        
 
"""
timestamps de la fenêtre de simulation + expressions compréhensibles pour humains
"""
tss= start + simStart * step
tse= start + simEnd * step
_s=tsToHuman(tss)
_e=tsToHuman(tse)


for i in range (int(simStart),int(simEnd)):
   RSB.SystemLoop(i,qdro_u,dromo)

integralite=False # booléen donnant l'intégralité ou non du graphique
if integralite:
    k1=simStart
    k2=simEnd
Tinj_dro=RSB.Tinj_dro[k1:k2]
Tsor_dro=RSB.Tsor_dro[k1:k2]
Tinj_sto=RSB.Tmoy_inj_sto[k1:k2]
Tsor_sto=RSB.Tmoy_sor_sto[k1:k2]
Tinj_pac=RSB.Tmoy_inj_pac[k1:k2]
Tsor_sto=RSB.Tmoy_sor_pac[k1:k2]

T1=Tsor_dro-Tinj_dro
T2=Tsor_sto-Tinj_sto
graphe(k1,k2)    
    

