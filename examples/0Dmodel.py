from dromosense import getCsvDatas
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

"""
col 0 : Température air (°C)
col 1 : Température point de rosée (°C)
col 2 : Nature des précipitations
col 3 : Vitesse du vent (m/s)
col 4 : Rayonnement global (W/m2)
col 5 : Rayonnement atmosphérique (W/m2)
"""
step, datas = getCsvDatas("meteo.csv",preview=True)


"""
albedo et émissivité
grandeurs sans unité
"""
albedo=0.08
epsilon=0.92

"""
L : Largeur de la chaussée en m
"""
L=4

"""
Hv Coefficient convectif entre la surface et l'air (fonction affine de la vitesse du vent)
en W/(m2K)
"""
Hv = 5.8 + 4.1*datas[:,3]
B1 = (1-albedo)*datas[:,4] + datas[:,5] + Hv*datas[:,0]

"""
conductivités thermiques des couches de surface, drainante et de base
unité : W/(m.K)
"""
ks=2.34
kd=1.56
kb=1.76

"""
capacités calorifiques volumiques des différentes couches
unité J/(m^3.K)
"""
Cs=2144309
Cd=1769723
Cb=2676728
Cf=4181000

"""
épaisseurs en m
"""
hs=0.06
hd= 0.08
hb=10

# constante de Stefan-Boltzmann en W/(m2K4)
#sigma=5.67*10**(-8)
sigma=5.67e-8

# température d'injection du fluide dans le dromotherme en °C
Tinjection=10

"""
on part du principe que l'on travaille à débit fixe de 50 l/h
"""
qf=0.05/3600 # en m^3/s

"""
surface d'échange entre le fluide et une paroi solide (couche drainante)
On suppose qu'elle est 100 fois plus grande de la surface du dromotherme
"""
S=L*100
"""
coefficient d'échange convectif entre le fluide et la couche drainante en W/(m2K)
perso, je trouve que c'est beaucoup, ce n'est quant même pas de l'eau bouillante ?
"""
h=2000

def rd(k1,k2,h1,h2):
    """
    calcule le coefficient d'échange surfacique entre 2 couches de conductivité k1 et k2 et d'épaisseurs h1 et h2 
    """
    return 2*k1*k2/(h1*k2+h2*k1)

"""
il faudrait les calculer dynamiquement à partir des valeurs de conductivité
"""
rds=27 # coefficient d'échange surfacique entre la couche drainante et la surface

rdb=3.28 # coefficient d'échange surfacique entre la couche drainante et la surface

print(rd(ks,kd,hs,hd))
print(rd(kd,kb,hd,hb))
input("press any key")

def Tf_out(Td):
    """
    calcule la température sortante du fluide pour une température de la couche drainante donnée
    varie selon la température d'injection....
    Td : scalaire ou vecteur
    """
    return Td+(Tinjection-Td)*math.exp(-h*S/(qf*Cf))


def F(t, X):
    """
    Définition du système d'équations différentielles 
    i is the time
    """
    Ts=X[0]
    Td=X[1]
    Tb=X[2]
    i=int(t/step)

    Tsortie=Tf_out(Td)
    
    y1 = ( B1[i] - Hv[i]*Ts - epsilon*sigma*(Ts+273.15)**4 - rds*(Ts-Td) ) / ( Cs*hs )
    y2 = ( rds*(Ts-Td) - rdb*(Td-Tb) - qf*Cf/L*(Tsortie-Tinjection) ) / ( Cd*hd )
    y3 = rdb*(Td-Tb) / ( Cb * hb )
    
    return [y1,y2,y3]


def G(X,t):
    """
    pour odeint fonction de l'ancienne API
    """
    return F(t,X)

"""
on recrée le temps discrétisé
"""
t=step*np.arange(datas.shape[0])

"""
cf https://docs.scipy.org/doc/scipy/reference/integrate.html
"""
solution = solve_ivp(F,[0,(datas.shape[0]-1)*step],[10,10,10],t_eval=t)
oldsol = odeint(G,[10,10,10],t)

Tsortie=Tf_out(solution.y[1])

figure1 = plt.figure(figsize = (10, 5))
plt.subplot(211)
plt.xlabel('Temps (en secondes)')
plt.ylabel('Températures (en °C)')
plt.title("Profils de températures 0D")
plt.plot(solution.t,solution.y[0],label="Température couche de surface")
plt.plot(solution.t,solution.y[1],label="Température couche drainante (solve_ivp new API)")
plt.plot(solution.t,solution.y[2],label="Température couche de base")
plt.plot(solution.t,Tsortie)
plt.legend(loc="upper right")
plt.subplot(212)
plt.plot(solution.t,solution.y[0],label="Température couche de surface")
plt.plot(solution.t,oldsol[:,1],label="Température couche drainante (odeint old API)")
plt.legend(loc="upper right")
plt.show()