import numpy as np
import matplotlib.pyplot as plt
from dromosense.tools import *
"""
permet d'importer automatiquement albedo, epsilon et sigma, kelvin, Cpf en J/(kg.K)
"""
from dromosense.constantes import *

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
meteo = np.loadtxt('meteo2.txt')
#meteo[:,1:]=meteo[0,1:]
f2 = 1000.0*1.1*(0.0036*meteo[:,4]+0.00423)
f1 = (1.0-albedo)*meteo[:,5] + meteo[:,6] + f2*(meteo[:,1]+kelvin)
t = meteo[:,0]
nt = t.size # le nombre de points dans la discretisation temporelle
dt = 3600.0

"""
l'initialisation complète de l'objet se fait de la manière suivante :
dromo=OneDModel('input.txt',dt,nt,L,dx,qf)
avec :
- dt : pas de temps en secondes
- nt : nombre de points dans la discrétisation temporelle
- L : largeur de la chaussée en m
- dx : pas d'espace en m
- qf : débit en m3
le fluide est de l'eau
la classe OneDModel va chercher les valeurs de Cpf et rho_eau dans le fichier des constantes
içi on ne donne pas L, dx ni qf car la classe dispose de valeurs par défaut (4,0.75,0/035/3600)
mais dans la pratique, il faut fournir les bonnes valeurs

ensuite :
- il faut affecter les vecteurs dromo.f1 et dromo.f2 issus de l'exploitation du monitoring :-)
- il faut fixer la condition initiale pour dromo.T[0,:,:] - cf plus loin en utilisant T2d[0,1:]

Nota : les variables de classe sont publiques en python, contrairement au C++
on y accède via nom_instance.nom_variable
içi l'instance est dromo
"""
# instanciation.....
dromo=OneDModel('input.txt',dt,nt)
#dromo.showVectors()
dromo.f1=f1
dromo.f2=f2
input("press any key")

# IMPORTATION DES DONNEES DES COUCHES DE CHAUSSEE
# on ne peut pas utiliser input comme nom de variable car c'est une fonction python
# ex : input("press any key")
_input = np.loadtxt('input.txt')
nc = _input.shape[0] # nombre de couches
ha = _input[:,0] # hauteur des couches en m
le = _input[:,1] # coef d'echanges des couches (derniere valeur non utilisee) en W/(m2.K)
rc = _input[:,2] # capacites calorifiques des couches en J/(m3.K)


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
labels=["T surface","T drainant", "T base 1", "T base 2", "T massif"]
T2d = kelvin + np.loadtxt('T2d2.txt')
plt.subplot(111)
plt.title("données du modèle 2D")
plt.ylabel("Température K")
plt.xlabel("temps en heures")
for i in range(1,T2d.shape[1]):
    plt.plot(T2d[:,i],label=labels[i-1])
plt.legend()
plt.show()



L = 4.0 # Largeur de chaussee en m
dx = 0.75 # pas d'espace en m
qf = 0.035/3600.0         # debit volumique du fluide (m^3/s)
Cf = Cpf*rho_eau #4200000.0 # capacite calorifique volumique de l'eau (J/(m3.K))
phi = 0.2     #( porosite de la couche drainante)
l=3 # largeur de la chaussée
nx = int(L/dx) # Nombre de blocs



"""
axe 0 : Temps
axe 1 : nombre de couches ou z
axe 2 : nombre de blocs ou x
"""
T = np.zeros((nt,nc,nx)) # Champ de temperature

# Conditions initiales
for i in range(nx) :
    T[0,:,i] = T2d[0,1:]
    dromo.T[0,:,i] = T2d[0,1:]

# Conditions aux limites
Tinj = 10.0 + kelvin
lambd = 10000000 # raideur pour imposer Tinj

# RESOLUTION DU CHAMP THERMIQUE
# Calcul des differents elements des vecteurs diagonaux (A), sur-diagonaux (B) et sous-diagonaux (C)
A = np.zeros((nc))
B = np.zeros((nc))
C = np.zeros((nc))

A[1:nc-1] = dt * (le[0:nc-2] + le[1:nc-1])
A[nc-1] = dt * le[nc-2]
#A[1] = A[1] + dt * (qf * Cf) / dx
A = A + ha*rc
print(A)
B[0:nc-1] = - dt * le[0:nc-1]
C[1:nc] = - dt * le[0:nc-1]

for n in range(1,nt):
    dromo.iterate(n,Tinj)
    for j in range(0,nx):
        A[0] = dt * (f2[n] + le[0] + 4.0*epsilon*sigma*T[n-1,0,j]**3) + ha[0] * rc[0]
        R = ha*rc*T[n-1,:,j]
        R[0] = R[0] + dt * (f1[n] + 3.0*epsilon*sigma*T[n-1,0,j]**4)
        if j==0:
           #R[1] = R[1] - dt * lambd * Tinj
           R[1] = dt*lambd*Tinj
           A[1] = dt*lambd*1.0
           C[1] = 0.0
           B[1] = 0.0
        else:
           R[1] = R[1] + dt * (qf * Cf / (l*dx)) * (T[n-1,1,j-1]-T[n-1,1,j])
           C[1] = - dt * le[0]
           B[1] = - dt * le[1]
           A[1] = dt * (le[0] + le[1]) + ha[1] * rc[1] #+ dt * (qf * Cf) / dx
        T[n,:,j] = sol_tridiag(A,B,C,R)

np.savetxt('T1d_summer.txt', T[:,:,-1]-kelvin, fmt='%.2e')

plt.subplot(111)
plt.title("modèle 1D vs modèle 2D")
plt.ylabel("°C")
plt.xlabel("heures")
plt.plot(t/3600,T[:,1,-1]-kelvin,label="1D model")
plt.plot(t/3600,T2d[:,2]-kelvin,'--',label="2D model")
plt.legend(loc="upper right")
plt.show()

plt.subplot(111)
plt.title("modèle 1D vs classe")
plt.ylabel("°C")
plt.xlabel("heures")
plt.plot(t/3600,T[:,1,-1]-kelvin,label="1D model")
plt.plot(t/3600,dromo.T[:,1,-1]-kelvin,label="classe")
plt.legend(loc="upper right")
plt.show()

import matplotlib as mpl

def colorFader(c1,c2,mix=0):
    """
    gradient entre les couleurs c1 et c2

    mix=0 : c1

    mix=1 : c2
    """
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

c1='#ffdede' #rouge très clair
c2='red'

plt.subplot(111)
plt.title("Modèle 1D - couche drainante \n gradient de température dans le fluide fonction de x\n profil en travers découpé en {} couches".format(nx))
plt.ylabel("°C")
plt.xlabel("heures")
for i in range(nx-1,-1,-1):
    plt.fill_between(t/3600,T[:,1,i]-kelvin,label="1D model",color=colorFader(c1,c2,i/nx))
plt.show()

plt.subplot(111)
plt.title("Modèle 1D \n Sortie de l'échangeur \n gradient de température fonction de z")
plt.ylabel("°C")
plt.xlabel("heures")
for i in range(0,nc):
    plt.fill_between(t/3600,T[:,i,-1]-kelvin,label="1D model",color=colorFader(c2,c1,i/nx))
plt.show()
