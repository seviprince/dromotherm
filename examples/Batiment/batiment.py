import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import cmath as math
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
# tu peux directement charger le fichier déposé par Frédéric, pas la peine de le recopier dans le répertoire courant....
# ../ veut dire que tu retournes un répertoire en arrière
# donc içi avec ../../ tu retournes au niveau de la racine de ton dossier github
meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')
#meteo = np.loadtxt('meteo3.txt')
# mieux vaut appeler cette variable Tconsigne, Tint pourrait servir pour héberger la température intérieure modélisée 
Tconsigne=19
#Tint=19
# tu n'es pas obligé de créer un vecteur spécifique...celà mobilise de la mémoire en plus
#Text=meteo[:,1]
#t=meteo[:,0]
Rm=8.24E-02 # Résistance thermique des murs (K/W)
Ri=1.43E-03 # Résistance superficielle intérieure
Rf=0.034 # Résistance due aux ifilitration et au renouvellement d'air
Rthe=1/(1/(Rm+Ri)+1/Rf) # Résistance thermique équivalente 
Ci=18.407 # Capacité thermique de l'air (Wh/K)
Cm=2636 # Capacité thermique des murs 
Besoin=(Tconsigne-meteo[:,1])/Rthe
#Graphique
figure1 = plt.figure(figsize = (10, 5))
ax1=plt.subplot(111)
plt.xlabel('Temps (en heure)')
plt.ylabel('Puissance (W)')
plt.title("Evolution des besoins d\'un bâtiment de 20m2 sur une année")
plt.plot(meteo[:,0],Besoin,label="besoin en Kwh",color="orange")
plt.legend(loc="upper left")
ax2=ax1.twinx()
plt.plot(meteo[:,0],meteo[:,1],label="Text")
plt.legend(loc="upper right")
plt.show()