import numpy as np
import matplotlib.pyplot as plt
from dromosense.tools import besoin_bat
from datetime import datetime
import time
from dateutil import tz
CET=tz.gettz('Europe/Paris')

"""
IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
0 : temps exprime en heure 
1 : temperature d'air (en deg Celsius)
2 : rayonnement global (en W/m2)
3 : rayonnement atmospherique (en W/m2)
4 : vitesse du vent (en m/s)
"""

def tsToTuple(ts):
    """
    ts : unix time stamp en s
    
    return date tuple tm_year, tm_mon, tm_mday, tm_hour, tm_min, tm_sec, tm_wday, tm_yday, tm_isdst
    """
    _time=datetime.fromtimestamp(ts,CET)
    _tuple=_time.timetuple()
    return(_tuple)


meteo = np.loadtxt('../../datas/corr1_RT2012_H1c_toute_annee.txt')
"""
cf https://www.unixtimestamp.com
01/01/2017 00:00:00 UTC+1 ou 01/01/2017 01:00:00 UTC
l'année est une hypothèse : Frédéric a téléchargé ces données sur le site de la RT2012 où la France est découpée en secteur
il s'agit de données fictives, compilées à partir à partir d'une dizaine d'années de données météo....
à tout à l'heure, peut-être ce sera plutôt 16h30 !
"""
start = 1483232400
summerStart = 1496278800
summerEnd = 1506819600

"""
building an agenda indicating wether people are working or not
1: work
0: rest
here we go for fixed working hours each day
start at 8 and stop at 17
"""
agenda=np.zeros(datas.shape[0])
time=start
tpl=tsToTuple(time)
work=0
if tpl.tm_hour in range(8,17):
    if tpl.tm_wday not in [5,6]:
        work=1
agenda[0]=work
for i in range (0,datas.shape[0]-1):
    tpl=tsToTuple(time)
    if tpl.tm_hour==17 and previous.tm_hour==16:
        if tpl.tm_wday not in [5,6]:
            work=0
    if tpl.tm_hour==8 and previous.tm_hour==7:
        if tpl.tm_wday not in [5,6]:
            work=1
    agenda[i]=work
    previous=tpl
    time+=step


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