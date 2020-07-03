import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from dromosense.tools import *
from dromosense.constantes import rho_eau,Cpf,kelvin
from scipy.integrate import odeint
#from scipy.integrate import solve_ivp
import math
from datetime import datetime
from dateutil import tz
CET=tz.gettz('Europe/Paris')

def ts_to_h(ts):
    """
    format a timestamp to something readable by a human
    """
    return datetime.fromtimestamp(ts,CET).strftime('%Y-%m-%d %H:%M:%S')

def StockLoop(i):
    """
    réalise une itération sur la température du stockage
    
    4 cas distincts :
    
    1) appel d'énergie en provenance du bâtiment + dromotherme en marche
    
    2) appel d'énergie en provenance du bâtiment + dromotherme arrêté
    
    3) pas d'appel d'énergie en provenance du bâtiment + dromotherme en marche
    
    4) pas d'appel d'énergie en provenance du bâtiment + dromotherme arrêté
    """
    pac=agenda_pac[i]
    dro=agenda_dro[i]
    if pac==1 and dro==1:
        der = (msto * cpsto * (Tinj_sto[i] - Tsor_sto[i]) - Pgeo[i]) / (m_sable * Cp_sable)
    if pac==1 and dro==0:
        der = - Pgeo[i] / (m_sable * Cp_sable)
    if pac==0 and dro==1:
        der = msto * cpsto * (Tinj_sto[i] - Tsor_sto[i]) / (m_sable * Cp_sable)
    if pac==0 and dro==0:
        der = 0
        
    diff[i+1]=der
    
    ## schéma de discrétisation
    return Tsable[i]+step*der

def SystemLoop(i):
    """
    On commence par mettre à jour les température d'injection et de sortie de la PAC
    
    On réalise ensuite une itération de dromotherme selon 2 cas distincts :
    
    - le dromotherme est en marche et le fluide circule avec un débit unitaire qdro_u
    
      on récupère de l'énergie et on alimente le stockage via l'échangeur de séparation de réseaux
      
    - le dromotherme est à l'arrêt : le débit est nul et l'échangeur de séparation de réseau ne tourne pas
    
      1) pas de prélèvement par l'échangeur de séparation de réseau
      
         la température d'injection dans le dromotherme est égale à la température de sortie
         
      2) fonctionnement à perte nulle
      
         les températures d'injection et de sortie au niveau du stockage sont égales à celles correspondant à l'itération précédante
         
    Dans tous les cas, on applique StockLoop
    """
    
    dro=agenda_dro[i]
    pac=agenda_pac[i]
    y = Tsable[i-1]
    if pac == 1 :
        Tinj_pac[i] = y-C*Pgeo[i]/k
        Tsor_pac[i] = Tinj_pac[i]-Pgeo[i]/(mpac*cpac)
    else:
        Tinj_pac[i] = Tinj_pac[i-1]
        Tsor_pac[i] = Tsor_pac[i-1]
    
    if dro == 1:
        dromo.iterate(i,Tinj_dro[i-1]+kelvin,qdro_u)
        Tsor_dro[i]=dromo.T[i,1,-1]-kelvin
        Tsor_sto[i] = ( k * y + B * Tsor_dro[i] ) / ( k + B)
        Tinj_sto[i] = Tsor_sto[i] + coeff * eff * (Tsor_dro[i] - Tsor_sto[i])
        Tinj_dro[i] = Tsor_dro[i] - eff * (Tsor_dro[i] - Tsor_sto[i])
        
    else:
        dromo.iterate(i,Tinj_dro[i-1]+kelvin,0)
        Tsor_dro[i]=dromo.T[i,1,-1]-kelvin
        Tinj_dro[i]=Tsor_dro[i]
        Tinj_sto[i] = Tinj_sto[i-1] 
        Tsor_sto[i] = Tsor_sto[i-1]
    
    Tsable[i]=StockLoop(i-1)

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
    if label:
        l1 = "{} / {}".format(l1,label)
    l2 = "température de consigne dans le bâtiment : {} °C".format(Tconsigne)
    plt.title("{}\n{} -> {}\n{}\n".format(l1,_s,_e,l2))
    
    ## graphe 1 - la route
    plt.ylabel('dromotherme')
    
    ax1.plot(xrange, agenda_dro[s:e], color=clearblue, label="dro ON/OFF")
    ax1.legend(loc="upper left")
    
    ax11 = ax1.twinx()
    ax11.plot(xrange, Tsor_dro[s:e], label="Tsor_dro", color="red")
    ax11.plot(xrange, Tinj_dro[s:e], label="Tinj_dro", color="purple")
    ax11.legend(loc="upper right")
    
    ## graphe 2 - la météo
    ax2 = plt.subplot(512, sharex=ax1)
    
    ax2.plot(xrange,meteo[s:e,2],label="rayonnement global en W/m2",color="orange")
    ax2.legend(loc="upper left")
    
    ax21 = ax2.twinx()
    ax21.plot(xrange, meteo[s:e,1], label="T ext")
    ax21.legend(loc="upper right")
    
    ## graphe 3 - le stockage
    ax3 = plt.subplot(513, sharex=ax1)
    plt.ylabel('stockage')
    
    ax3.plot(xrange, diff[s:e], color=clearblue, label="derivée")
    ax3.legend(loc="upper left")
    
    ax31 = ax3.twinx()
    ax31.plot(xrange, Tinj_sto[s:e], label="Tinj_sto", color="orange")
    ax31.plot(xrange, Tsor_sto[s:e], label="Tsor_sto", color="blue")
    ax31.plot(xrange, Tsable[s:e], label="Tsable", color="red")
    ax31.legend(loc="upper right")
    
    ## graphe 4 - la PAC
    ax4 = plt.subplot(514, sharex=ax1)
    plt.ylabel('PAC')
    
    ax4.plot(xrange, Tinj_pac[s:e], label="Tinj_pac", color="red")
    ax4.plot(xrange, Tsor_pac[s:e], label="Tsor_pac", color="#7cb0ff")
    ax4.plot(xrange, Tinj_pac[s:e]-Tsor_pac[s:e], label="deltaT PAC", color="k")
    ax4.legend()
    
    ## graphe 5 - le besoin du bâtiment
    ax5 = plt.subplot(515, sharex=ax1)
    plt.ylabel('Bâtiment')
    plt.xlabel("Temps - 1 unité = {} s".format(step))
    
    ax5.plot(xrange, agenda_pac[s:e], color=clearblue, label="PAC ON/OFF")
    ax5.legend(loc="upper left")
    
    ax51 = ax5.twinx()
    ax51.plot(xrange, besoin_total[s:e], label="besoin {} en W".format(heating), color="orange")
    ax51.legend(loc="upper right")
    
    plt.show()

def ECSPower(min, max):
    """
    on modélise l'eau du réseau comme une fonction sinusoidale de période annuelle
    
    cette fonction est complètement calée sur un fichier météo qui commence au 1er janvier mais qui peut être pluriannuel
    
    min : température minimale d'injection de l'eau du réseau dans le ballon

    max : température maximale d'injection de l'eau du réseau dans le ballon
    """
    T_water=np.zeros(meteo.shape[0])
    ## période
    w = 2*math.pi/npy
    for i in range(meteo.shape[0]):
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
"""

#meteo = np.loadtxt('../../datas/meteo_Bourget_du_lac.txt') 
meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')
npy = meteo.shape[0]
meteo=np.concatenate((meteo,meteo))
meteo[npy:meteo.shape[0],0]=meteo[npy:meteo.shape[0],0]+npy
print(meteo.shape)
f2 = 1000.0*1.1*(0.0036*meteo[:,4]+0.00423)
f1 = (1.0-albedo)*meteo[:,2] + meteo[:,3] + f2*(meteo[:,1]+kelvin)

# longueur et largeur de l'échangeur, exprimées en m
lincha=7.5
larcha=4
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

# température d'entrée et de sortie du fluide dans le stockage
Tinj_sto=np.zeros(meteo.shape[0])
# température de sortie du fluide après transit dans le stockage
Tsor_sto=np.zeros(meteo.shape[0])
# température d'entrée du fluide géothermique dans le stockage (sortie de la PAC)
Tsor_pac=np.zeros(meteo.shape[0])
# température de sortie du fluide géothermique dans le stockage (entrée  de la PAC)
Tinj_pac=np.zeros(meteo.shape[0])

# Tsable : Température du stockage/sable
Tsable=np.zeros(meteo.shape[0])
# valeur de la dérivée de la température du stockage
diff=np.zeros(meteo.shape[0])

"""
On initialise les températures d'injection et de sortie dans le dromotherme à 10
de fait, quant il n'y aura pas de récupération d'énergie, le bilan énergétique sera nul
"""
# température d'injection du fluide dans le dromotherme
Tinj_dro=10*np.ones(meteo.shape[0])
# température de sortie de fluide après collecte de chaleur dans le dromotherme
Tsor_dro=10*np.ones(meteo.shape[0])


"""
massif de stockage
"""
m_sable=70200.0 # masse de sable en kg
Cp_sable=1470.0 # capacité calorifique massique du sable en J/Kg.K

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

print("le k du système géothermique vaut {} W/K".format(k))

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

B = (msto * cpsto -k/2) * coeff * eff

print("coeff vaut {} B vaut {} W/K".format(coeff, B))

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
C=1-k/(2*mpac*cpac)

"""
on calcule les index permettant de boucler sur l'été, ainsi que le besoin net du bâtiment et la puissance géothermique à développer
"""
i_summerStart=(summerStart-start)//step

i_summerEnd=i_summerStart+(summerEnd-summerStart)//step

apport_solaire = Scap * FSm * meteo[:,2]

besoinBrut = besoin_bat(Tconsigne,meteo[:,1],Rm,Ri,Rf)

besoin_chauffage = besoinBrut - apport_solaire

"""
on efface les besoins de chauffage sur les étés
"""

besoin_chauffage[i_summerStart:i_summerEnd] = np.zeros(i_summerEnd-i_summerStart)

besoin_chauffage[i_summerStart+npy:i_summerEnd+npy]=np.zeros(i_summerEnd-i_summerStart)

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

ECS = ECSPower(Tentree_hiver, Tentree_ete)

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
ECSupply=False
usecase=3
"""
*************************************
*************************************
"""
            
if ECSupply:
    heating = "chauffage+ECS"
else:
    heating = "chauffage"
            
besoin_total = besoin_chauffage + ECSupply * ECS

besoin_surfacique = besoin_total / Scap

Pgeo = (COP-1) * besoin_total / COP

label = ""


# initialisation des agendas à 0 : aucun équipement en fonctionnement par défaut
agenda_dro=np.zeros(meteo.shape[0])
agenda_pac=np.zeros(meteo.shape[0])

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
    simEnd = i_summerEnd
    agenda_dro[simStart:simEnd]=np.ones(simEnd-simStart)
    agenda_pac[simStart:simEnd]=np.ones(simEnd-simStart)
    heating = "ECS"

if usecase == 1:
    # dromotherme durant l'été
    simStart = i_summerStart
    simEnd=i_summerEnd+2000
    agenda_dro[simStart:i_summerEnd]=np.ones(i_summerEnd-simStart)
    agenda_pac[simStart:simEnd]=np.ones(simEnd-simStart)

if usecase == 2:
    # simulation annuelle
    simStart = i_summerStart
    simEnd=i_summerStart+365*24
    agenda_dro[simStart:simEnd]=np.ones(simEnd-simStart)
    agenda_pac[simStart:simEnd]=np.ones(simEnd-simStart)

if usecase == 3:
    # simulation annuelle
    # dromotherme l'été et par intermittence l'hiver quant le rayonnement global est au dessus de 250 W/m2
    simStart = i_summerStart
    simEnd=i_summerStart+365*24
    agenda_dro[simStart:i_summerEnd]=np.ones(i_summerEnd-simStart)
    agenda_pac[simStart:simEnd]=np.ones(simEnd-simStart)
    level=250
    label = "hiver si ray.>250 W/m2 : dromo on".format(level) 
    for i in range(simStart,simEnd):
        if meteo[i,2] >= level:
            agenda_dro[i]=1


for i in range(simStart,simEnd):
    if Pgeo[i]==0:
        agenda_pac[i]=0            

"""
timestamps de la fenêtre de simulation + expressions compréhensibles pour humains
"""
tss= start + simStart * step
tse= start + simEnd * step
_s=ts_to_h(tss)
_e=ts_to_h(tse)


for i in range (int(simStart),int(simEnd)):
   SystemLoop(i)

"""
BILAN ENERGETIQUE
"""
Surface_dro=larcha*lincha
# d et f : index de début et de fin sur lesquels on va réaliser le bilan
# dr et fr : index de début et de fin de la récupération/collecte énergétique
d = i_summerStart
f = simEnd
dr = i_summerStart
fr = i_summerEnd
Pdro=mdro*cpdro*(Tsor_dro[d:f]-Tinj_dro[d:f])/Surface_dro # en W/m^2
# toutes valeurs en kWh/m^2
Edro=(np.sum(Pdro))*step/(3600*1000)
Esolaire=(np.sum(meteo[dr:fr,2]))*step/(3600*1000)
conso_bat=(np.sum(besoin_surfacique[d:f]))*step/(3600*1000)
Eelec=conso_bat/COP
Eprimaire=2.58*Eelec
# taux de récupération en %
Taux=Edro*100/Esolaire

"""
affichages énergétiques
"""
print("Bilan énergétique sur la période : {} à {}".format(_s,_e))
print('Energie récupérée par le dromotherm : {} kWh/m2'.format(Edro))
print('Energie solaire recue : {} kWh/m2'.format(Esolaire))
print('Taux de récupération : {} %'.format(Taux))
print('Consommation du bâtiment : {} kWh/m2'.format(conso_bat))
print('Energie électrique consommée par la PAC : {} kWh/m2'.format(Eelec))
print('Energie primaire consommée par la PAC : {} kWh/m2'.format(Eprimaire))
    
graphe(simStart,simEnd)