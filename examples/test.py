from dromosense.constantes import *
import numpy as np

print(albedo)
print(Cs)

schedule=np.array([[8,17],[8,17],[8,17],[8,17],[8,17],[-1,-1],[-1,-1]])

meteo = np.loadtxt('../datas/corr1_RT2012_H1c_toute_annee.txt')

start = 1483232400
summerStart = 1496278800
#summerEnd = 1506819600 : On s'était trompé sur la période de récupération de chaleur que Monsieur Frédéric a envoyée..c'est du 1er juin au 31 août et non au 30 septembre
summerEnd=1504141200
step=3600

nbpts=meteo.shape[0]

from dromosense.tools import *

agenda=basicAgenda(nbpts,step,start,summerStart,summerEnd,schedule=schedule)

print(agenda.shape)
print(meteo.shape)

plt.subplot(111)
plt.plot(agenda)
plt.show()