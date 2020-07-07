import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from dromosense.tools import *
from dromosense.constantes import rho_eau,Cpf,kelvin

start = 1483232400
summerStart = 1496278800
#summerEnd = 1506819600 # 30 septembre
summerEnd=1504141200 # 30 août
step=3600

# débit dans le dromotherme
qdro = 0.035/step # m3/s

"""
IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
0 : temps exprime en heure 
1 : temperature d'air (en deg Celsius)
2 : rayonnement global (en W/m2)
3 : rayonnement atmospherique (en W/m2)
4 : vitesse du vent (en m/s)
"""
meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')
print(meteo.shape)
f2 = 1000.0*1.1*(0.0036*meteo[:,4]+0.00423)
f1 = (1.0-albedo)*meteo[:,2] + meteo[:,3] + f2*(meteo[:,1]+kelvin)

dromo=OneDModel('input.txt',step,meteo.shape[0],4,0.75)
dromo.f1 = f1
dromo.f2 = f2

Tinj=10+kelvin

i_summerStart=(summerStart-start)//step
i_summerEnd=i_summerStart+(summerEnd-summerStart)//step

dromo.T[i_summerStart,:,:] = np.ones((dromo.T.shape[1],dromo.T.shape[2]))*10+kelvin
print("activation du dromotherm a l'index {} avec la condition initiale suivante :".format(i_summerStart))
print(dromo.T[i_summerStart,:,:])
input("press any key")

for n in range(i_summerStart,i_summerEnd):
    dromo.iterate(n,Tinj,qdro)

mdro = qdro * rho_eau
cpdro = Cpf
Pdro = mdro * cpdro * np.maximum(dromo.T[:,1,-1] - Tinj , np.zeros(meteo.shape[0]))

Edro=np.sum(step*Pdro)/(3600*1000)
Esolaire=np.sum(step*4*meteo[:,2][i_summerStart:i_summerEnd])/(3600*1000)
print("Energie récupérée par le dromotherme {}".format(Edro))
print("Energie reçue par le dromotherme {}".format(Esolaire))
# taux de récupération en %
Taux=Edro*100/Esolaire
print("rendement {}".format(int(Taux)))

plt.subplot(211)
plt.plot(dromo.T[:,1,-1]-kelvin,label="1D model T sortie couche drainante °C")
plt.legend()
plt.subplot(212)
plt.plot(4*meteo[:,2], label="rayonnement sur dromotherme : {} W".format(int(Esolaire)), color="orange")
plt.plot(Pdro, label="énergie captée par le dromotherme : {} W soit un rendement de {} %".format(int(Edro),int(Taux)), color="green")
plt.legend()
plt.show()