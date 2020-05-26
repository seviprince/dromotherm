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
def besoin_bat(Tint,Text,Rm,Ri,Rf):
    Rthe=1/(1/(Rm+Ri)+1/Rf)# Résistance thermique équivalente 
    return (Tint-Text)/Rthe
meteo = np.loadtxt('meteo3.txt')
Tint=19
Text=meteo[:,1]
t=meteo[:,0]
R=np.array(3)
C=np.array(2)
Rm=8.24E-02 # Résistance thermique des murs (K/W)
Ri=1.43E-03 # Résistance superficielle intérieure
Rf=0.034 # Résistance due aux infiltrations+vitre et au renouvellement d'air 
Ci=18.407 # Capacité thermique de l'air (Wh/K)
Cm=2636 # Capacité thermique des murs 
Besoin= besoin_bat(Tint,Text,Rm,Ri,Rf)
#Graphique
figure1 = plt.figure(figsize = (10, 5))
plt.xlabel('Temps (en heure)')
plt.ylabel('Puissance (W)')
plt.title("Evolution des besoins d\'un bâtiment de 20m2 sur une année")
plt.plot(t,Besoin)
plt.show()