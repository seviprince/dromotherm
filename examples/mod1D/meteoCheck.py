import numpy as np
import matplotlib.pyplot as plt

"""
what is what in meteo files ?
"""

"""
IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
0 : temps exprime en heure 
1 : temperature d'air (en deg Celsius)
2 : rayonnement global (en W/m2)
3 : rayonnement atmospherique (en W/m2)
4 : vitesse du vent (en m/s)
"""
meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')

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
meteo2 = np.loadtxt('meteo2.txt')

labels=["T surface","T drainant", "T base 1", "T base 2", "T massif"]
T2d = np.loadtxt('T2d2.txt')
offset=3599
t=np.arange(offset,offset+meteo2.shape[0])
plt.subplot(111)
plt.plot(t,T2d[:,1],label="T drainant > T2d.txt")
plt.plot(meteo[:,1],label="Text > datas/corr1_RT2012_H1c_toute_annee.txt")
plt.plot(t,meteo2[:,1],'.',label="Text > meteo2.txt")
plt.legend()
plt.show()

