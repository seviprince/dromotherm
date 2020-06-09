import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from dromosense.tools import *
from dromosense.constantes import rho_eau,Cpf,kelvin
from scipy.integrate import odeint
#from scipy.integrate import solve_ivp
import cmath as math

# débit dans le dromotherme
qdro = 0.035/3600 # m3/s

start = 1483232400
summerStart = 1496278800
#summerEnd = 1506819600 # 30 septembre
summerEnd=1504141200 # 30 août
step=3600

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

dromo=OneDModel('input.txt',step,meteo.shape[0],4,0.75,qdro)
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
    dromo.iterate(n,Tinj)
    
"""
T1d contiendra les températures des couches à la sortie du dromotherme
pour charger T1d.txt
_test = np.loadtxt('T1d.txt')
axe 0 : temps
axe 1 : z
_test[:,0] : couche de surface
_test[:,1] : couche drainante
"""
np.savetxt('T1d.txt', dromo.T[:,:,-1]-kelvin, fmt='%.2e')
    
plt.subplot(111)
plt.plot(dromo.T[:,1,-1]-kelvin,label="1D model T sortie couche drainante")
plt.legend()
plt.show()