import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from dromosense.tools import *
from dromosense.constantes import rho_eau,Cpf,kelvin
from scipy.integrate import odeint
#from scipy.integrate import solve_ivp
import math

verbose = False

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
print(dromo.T[0,:,:])


def F(y,t):
    """
    y : Tsable = Tstockage
    
    t : time index
    
    result : dy/dt = dTsable/dt = dTstockage/dt
    """
    i = int(t) 
    
    if verbose:
        print("we have t={} and y={}".format(i,y))
    
    # Je mets directement la forme développée de dy/dt sans calculs intermédiaires qui risqueraient de "forcer" un couplage      
    der = (dro*msto * cpsto * (((( k * y + B * Tsor_dro[i] ) / ( k + B)) + coeff * eff * (Tsor_dro[i] - (( k * y + B * Tsor_dro[i] ) / ( k + B)))) -(( k * y + B * Tsor_dro[i] ) / ( k + B)) ) - Pgeo[i] * pac ) / (m_sable * Cp_sable)
               
    return der

# température d'entrée et de sortie du fluide dans le stockage
Tinj_sto=np.zeros(meteo.shape[0])
# température de sortie du fluide après transit dans le stockage
Tsor_sto=np.zeros(meteo.shape[0])
# température d'entrée du fluide géothermique dans le stockage (sortie de la PAC)
Tsor_pac=np.zeros(meteo.shape[0])
# température de sortie du fluide géothermique dans le stockage (entrée  de la PAC)
Tinj_pac=np.zeros(meteo.shape[0])
# température du stockage
Tsable=np.zeros(meteo.shape[0])


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

print("nous allons simuler la récolte énergétique entre les heures {} et {}".format(i_summerStart,i_summerEnd))

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
"""
on modélise l'eau du réseau comme une fonction sinusoidale de période annuelle
"""
w=2*math.pi/npy
T_eau=np.zeros(meteo.shape[0])
besoin_ECS=np.zeros(meteo.shape[0])
for i in range(i_summerStart,i_summerStart+npy):
    #T_eau[i]=(1 + math.sin(w*(i-summerStart))) * (Tentree_ete-Tentree_hiver) / 2
    T_eau[i]= ((Tentree_ete-Tentree_hiver) / 2)* math.sin(w*(i-summerStart)) + (Tentree_ete+Tentree_hiver) / 2
    # le besoin s'entend pour une journée, ie 24*3600 secondes
    # il faut donc diviser par 24*3600 pour convertir de J à W, Cpf étant exprimée en J/kg/K
    besoin_ECS[i]=Volume_ballon*Npers*(Tballon-T_eau[i])*Cpf/(24*3600) 
    

besoin_total=besoin_chauffage+besoin_ECS

besoin_surfacique=besoin_total/Scap

Pgeo=(COP-1)*besoin_total/COP


"""
SOLVEUR
Tsable : Température du stockage/sable
"""
# changer usecase pour tester différentes choses
usecase=3

from datetime import datetime
from dateutil import tz
CET=tz.gettz('Europe/Paris')

def ts_to_h(ts):
    return datetime.fromtimestamp(ts,CET).strftime('%Y-%m-%d %H:%M:%S')


_s=ts_to_h(summerStart)


# initialisation des agendas à 0 : aucun équipement en fonctionnement par défaut
agenda_dro=np.zeros(meteo.shape[0])
agenda_pac=np.zeros(meteo.shape[0])

if usecase == 1:
    simEnd=i_summerEnd+4000
    _e=ts_to_h(summerStart+(simEnd-i_summerStart)*step)
    # dromotherme durant l'été et chauffage à partir du stock en continu durant 2000 heures en suivant
    label="dromotherme durant l'été et chauffage à partir du stock jusqu'au {}".format(_e)
    agenda_dro[i_summerStart:i_summerEnd]=np.ones(i_summerEnd-i_summerStart)
    agenda_pac[i_summerEnd:simEnd]=np.ones(simEnd-i_summerEnd)

if usecase == 2:
    simEnd=i_summerStart+365*24
    _e=ts_to_h(summerStart+(simEnd-i_summerStart)*step)
    # simulation annuelle
    label="dromotherme et PAC sur ON toute l'année "
    agenda_dro[i_summerStart:simEnd]=np.ones(simEnd-i_summerStart)
    agenda_pac[i_summerStart:simEnd]=np.ones(simEnd-i_summerStart)

if usecase == 3:
    simEnd=i_summerStart+365*24
    _e=ts_to_h(summerStart+(simEnd-i_summerStart)*step)
    # simulation à l'année
    # dromotherme l'été et par intermittence l'hiver quant le rayonnement global est au dessus de 250 W/m2
    agenda_dro[i_summerStart:i_summerEnd]=np.ones(i_summerEnd-i_summerStart)
    agenda_pac[i_summerEnd:simEnd]=np.ones(simEnd-i_summerEnd)
    for i in range(i_summerEnd,simEnd):
        if meteo[i,2] >= 250:
            label="dromotherme l'été et l'hiver quant le rayonnement global est au dessus de 250 W/m2"
            agenda_dro[i]=1
            


input("press any key")
plt.subplot(211)
plt.plot(agenda_dro,label="fonctionnement dromotherme")
plt.legend()
plt.subplot(212)
plt.plot(agenda_pac,label="fonctionnement pac")
plt.legend()
plt.show()

#Tsable = odeint(F,10,meteo[i_summerStart:simEnd,0]*step)

for i in range (int(i_summerStart),int(simEnd)):
    

    dro=agenda_dro[i]
    pac=agenda_pac[i]
      
      
    """
    SI dro==1
        - dromotherme en fonctionnement, on récupère de l'énergie et on alimente le stockage via l'échangeur de séparation de réseaux
        - si pac==1, on en tient compte dans le calcul de der en incluant une consommation à hauteur de Pgeo[i]
    SINON
        - dromotherme arrêté, pas de fonctionnement ni du dromotherme, ni de l'échangeur de séparation de réseau....
        - donc pas d'alimentation du stockage côté dromotherme
        - si pac==1, on en tient compte dans le calcul de der en incluant une consommation à hauteur de Pgeo[i]
    """
    
    if dro == 1:
        dromo.iterate(i,Tinj_dro[i-1]+kelvin,qdro_u)
        
        Tsor_dro[i]=dromo.T[i,1,-1]-kelvin
        
        Tsable[i]=Tsable[i-1]+step*F(Tsable[i-1],i) # La ligne clée du code: on utilise un Euler explicite pour déterminer Tsable ; un  Euler implicite serait un peu compliquée
    
        Tsor_sto[i] = ( k * Tsable[i]  + B * Tsor_dro[i] ) / ( k + B)
    
        Tinj_sto[i] = (( k * Tsable[i] + B * Tsor_dro[i] ) / ( k + B)) + coeff * eff * (Tsor_dro[i] - (( k * Tsable[i] + B * Tsor_dro[i] ) / ( k + B)))
    
        Tinj_dro[i] = Tsor_dro[i] - eff * (Tsor_dro[i] - Tsor_sto[i])
       
        
    else:
        
        dromo.iterate(i,Tinj_dro[i-1]+kelvin,0)
        
        Tsor_dro[i]=dromo.T[i,1,-1]-kelvin
        
        Tinj_dro[i]=Tsor_dro[i]

        Tinj_sto[i] = Tinj_sto[i-1] 
        
        Tsor_sto[i] = Tsor_sto[i-1]
             
    
    """s
    Si la PAC fonctionne, on met à jour les températures d'entrée et de sortie de PAC
    """ 
    
    if pac == 1 :
        
        Tinj_pac[i] = Tsable[i]-C*Pgeo[i]/k

        Tsor_pac[i] = Tinj_pac[i]-Pgeo[i]/(mpac*cpac)
    else:
       
        
        Tinj_pac[i] = Tinj_pac[i-1]
        Tsor_pac[i] = Tsor_pac[i-1]
    
        

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

"""
courbes et graphiques
"""
figure = plt.figure(figsize = (10, 10))
matplotlib.rc('font', size=8)

def graphe(start,stop):
    
    ax1 = plt.subplot(511)
    l1="couplage dromotherme/échangeur de séparation de réseau/stockage/PAC et simulation été/hiver"
    l2="température de consigne dans le bâtiment : {} °C".format(Tconsigne)
    plt.title("{}\n{}\n{}\n".format(l1,l2,label))
    plt.ylabel('dromotherme °C')
    ax1.plot(meteo[start:stop,0],Tsor_dro[start:stop],label="Tsor_dro",color="red")
    ax1.plot(meteo[start:stop,0],Tinj_dro[start:stop],label="Tinj_dro",color="purple")
    ax1.legend(loc="upper left")
    ax11=ax1.twinx()
    ax11.plot(agenda_dro,label="dro ON/OFF")
    ax11.legend(loc="upper right")
    
    ax2 = plt.subplot(512, sharex=ax1)
    ax2.plot(meteo[start:stop,0],meteo[start:stop,2],label="rayonnement global en W/m2",color="orange")
    ax22=ax2.twinx()
    ax22.plot(meteo[start:stop,0],meteo[start:stop,1],label="T ext")
    ax2.legend(loc="upper left")
    ax22.legend(loc="upper right")
    
    ax3 = plt.subplot(513, sharex=ax1)
    plt.ylabel('stockage °C')
    ax3.plot(meteo[start:stop,0],Tinj_sto[start:stop],label="Tinj_sto",color="orange")
    ax3.plot(meteo[start:stop,0],Tsor_sto[start:stop],label="Tsor_sto",color="blue")
    ax3.plot(meteo[start:stop,0],Tsable[start:stop],label="Tsable",color="red")
    ax3.legend(loc="upper left")
    #ax21=ax2.twinx()
   # ax21.plot(diff,label="derivée",color="black")
    #ax21.legend(loc="upper right")
    
    
    
    ax4 = plt.subplot(514, sharex=ax1)
    plt.ylabel('PAC °C')
    ax4.plot(meteo[start:stop,0],Tinj_pac[start:stop],label="Tinj_pac",color="red")
    ax4.plot(meteo[start:stop,0],Tsor_pac[start:stop],label="Tsor_pac",color="#7cb0ff")
    ax4.plot(meteo[start:stop,0],Tinj_pac[start:stop]-Tsor_pac[start:stop],label="ecart de températures de la PAC",color="k")
    ax4.legend()
    
    ax5 = plt.subplot(515, sharex=ax1)
    plt.ylabel('Bâtiment W')
    plt.xlabel("Temps - 1 unité = {} s".format(step))
    ax5.plot(meteo[start:stop,0],besoin_total[start:stop],label="besoin total du bâtiment net W",color="orange")
    ax5.plot(meteo[start:stop,0],besoin_ECS[start:stop],label="besoin en ECS ",color="g")
    ax5.plot(meteo[start:stop,0],besoin_chauffage[start:stop],label="besoin chauffage W",color="red")
    ax5.legend()
    
    plt.show()


integralite=True
if integralite:
    start=0
    stop=meteo.shape[0]
else:
    start=i_summerStart
    stop=simEnd
    
graphe(start,stop)






