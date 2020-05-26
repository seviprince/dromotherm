import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import cmath as math

"""
IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
0 : temps exprime en heure 
1 : temperature d'air (en deg Celsius)
2 : rayonnement global (en W/m2)
3 : rayonnement atmospherique (en W/m2)
4 : vitesse du vent (en m/s)
"""
<<<<<<< HEAD
def besoin_bat(Tint,Text,Rm,Ri,Rf):
    Rthe=1/(1/(Rm+Ri)+1/Rf)# Résistance thermique équivalente 
    return (Tint-Text)/Rthe
meteo = np.loadtxt('meteo3.txt')
Tint=19
Text=meteo[:,1]
t=meteo[:,0]
R=np.array(3)
C=np.array(2)
=======
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
>>>>>>> efe92150563d4d2bcae6b60dadadd9be2c758421
Rm=8.24E-02 # Résistance thermique des murs (K/W)
Ri=1.43E-03 # Résistance superficielle intérieure
Rf=0.034 # Résistance due aux infiltrations+vitre et au renouvellement d'air 
Ci=18.407 # Capacité thermique de l'air (Wh/K)
Cm=2636 # Capacité thermique des murs 
<<<<<<< HEAD
Besoin= besoin_bat(Tint,Text,Rm,Ri,Rf)
=======
Besoin=(Tconsigne-meteo[:,1])/Rthe
>>>>>>> efe92150563d4d2bcae6b60dadadd9be2c758421
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