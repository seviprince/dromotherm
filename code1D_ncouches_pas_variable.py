import numpy as np
import matplotlib.pyplot as plt

path="examples/mod1D"

#IMPORTATION DES DONNEES DES COUCHES DE CHAUSSEE
_input = np.loadtxt('{}/input.txt'.format(path))
nc = _input.shape[0] # nombre de couches
ha = _input[:,0] # hauteur des couches
le = _input[:,1] # coef d'echanges des couches (derniere valeur non utilisee)
rc = _input[:,2] # capacites calorifiques des couches

kelvin = 273.15

#IMPORTATION DES TEMPERATURES DU MODELE 2D
T2d = kelvin + np.loadtxt('{}/T2d2.txt'.format(path))


albedo = 0.08
epsilon = 0.92 
sigma = 5.67e-8

#IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
#1 : temps exprime en secondes depuis le 01/06 00:00
#2 : temperature d'air (en deg Celsius)
#3 : temperature du point de rosee (pas utile)
#4 : nature des precipitations (pas utile)
#5 : vitesse du vent (en m/s)
#6 : rayonnement global (en W/m2)
#7 : rayonnement atmospherique (en W/m2)
meteo = np.loadtxt('{}./meteo2.txt'.format(path))
f2 = 1000.0*1.1*(0.0036*meteo[:,4]+0.00423)   
f1 = (1.0-albedo)*meteo[:,5] + meteo[:,6] + f2*(meteo[:,1]+kelvin)
print(f1)
t = meteo[:,0]

# MAILLAGE X
#x = np.loadtxt('maillage_x.txt')
#x = np.array([0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.50,3.75,4])
x = np.array([0,0.75,1.5,2.25,3.0,3.75])
#nx = 10
nx = x.shape[0]-1
dx = 0.75*np.ones((nx))
#nx = x.shape[0] # nombre de points en x
#L = x[nx-1]
#dx = x[1:]-x[0:nx-1]
print(dx)

#L = 4.0 # Largeur de chaussee
#dx = 0.75 # pas d'espace
qf = 0.035/3600.0         # debit volumique du fluide (m^3/s)
Cf = 4200000.0 # capacite calorifique volumique de l'eau
phi = 0.0     #( porosite de la couche drainante)

dt = 3600.0

#nx = int(L/dx) # Nombre de blocs
nt = t.size # le nomnre de points dans la discretisation temporelle

"""
axe 0 : Temps
axe 1 : nombre de couches ou z
axe 2 : nombre de blocs ou x
"""
T = np.zeros((nt,nc,nx)) # Champ de temperature

# Conditions initiales
for i in range(nx) :
    T[0,:,i] = T2d[0,1:]

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


def sol_tridiag(A,B,C,D):
    N = A.size
    alpha=np.zeros((N))
    beta=np.zeros((N))
    X=np.zeros((N))
    alpha[0]=A[0]

    beta[0]=D[0]/alpha[0]
    for i in range(0,N-1):
        alpha[i+1]=A[i+1]-(C[i+1]*B[i])/alpha[i]
        beta[i+1]=(D[i+1]-B[i]*beta[i])/alpha[i+1]
    X[N-1]=beta[N-1]
    for i in range(N-2,-1,-1):
        X[i]=beta[i]-(B[i]*X[i+1]/alpha[i])
    return X


for n in range(1,nt):
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
           R[1] = dx[j-1]*R[1] + dt * (qf * Cf) * (T[n-1,1,j-1]-T[n-1,1,j])
           C[1] = - dx[j-1]*dt * le[0]
           B[1] = - dx[j-1]*dt * le[1]
           A[1] = dx[j-1]*dt * (le[0] + le[1]) + dx[j-1]*ha[1] * rc[1]
        T[n,:,j] = sol_tridiag(A,B,C,R)


np.savetxt('T1d.txt', T[:,:,-1]-kelvin, fmt='%.2e')

print(T[0,nx-1,:])
plt.plot(t,T[:,1,-1]-kelvin,label="1D model")
plt.plot(t,T2d[:,2]-kelvin,'r-',label="T2d")
plt.legend()
plt.show()
