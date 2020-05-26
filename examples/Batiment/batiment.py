import numpy as np
import matplotlib.pyplot as plt
from dromosense.tools import besoin_bat
"""
IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
0 : temps exprime en heure 
1 : temperature d'air (en deg Celsius)
2 : rayonnement global (en W/m2)
3 : rayonnement atmospherique (en W/m2)
4 : vitesse du vent (en m/s)
"""


meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')
Tconsigne=19
Rm=8.24E-02 # Résistance thermique des murs (K/W)
Ri=1.43E-03 # Résistance superficielle intérieure
Rf=0.034 # Résistance due aux infiltrations+vitre et au renouvellement d'air 
Ci=18.407 # Capacité thermique de l'air (Wh/K)
Cm=2636 # Capacité thermique des murs 
Besoin= besoin_bat(Tconsigne,meteo[:,1],Rm,Ri,Rf)

#Graphique
figure1 = plt.figure(figsize = (10, 5))
ax1=plt.subplot(111)
plt.xlabel('Temps (en heure)')
plt.ylabel('Puissance (W)')
plt.title("Evolution des besoins d\'un bâtiment de 20m2 sur une année")
plt.plot(meteo[:,0],Besoin,label="besoin en W",color="orange")
plt.legend(loc="upper left")
ax2=ax1.twinx()
plt.plot(meteo[:,0],meteo[:,1],label="Text en °C")
plt.legend(loc="upper right")
plt.show()