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
meteo = np.loadtxt('meteo3.txt')
Tint=19
Text=meteo[:,1]
t=meteo[:,0]
Rm=8.24E-02 # Résistance thermique des murs (K/W)
Ri=1.43E-03 # Résistance superficielle intérieure
Rf=0.034 # Résistance due aux ifilitration et au renouvellement d'air
Rthe=1/(1/(Rm+Ri)+1/Rf) # Résistance thermique équivalente 
Ci=18.407 # Capacité thermique de l'air (Wh/K)
Cm=2636 # Capacité thermique des murs 
Besoin=(Tint-Text)/Rthe
#Graphique
figure1 = plt.figure(figsize = (10, 5))
plt.xlabel('Temps (en heure)')
plt.ylabel('Puissance (W)')
plt.title("Evolution des besoins d\'un bâtiment de 20m2 sur une année")
plt.plot(t,Besoin)
plt.show()