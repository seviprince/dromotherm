import math
import numpy as np
import matplotlib.pyplot as plt

# pas de temps en secondes
step=3600
npy=24*365
# conductivité en W/(K.m)
lambda_sable=1.7
# masse volumique du sable en kg/m3
rho_sable=1700
# capacité calorifique massique du sable humide en J/Kg.K
cpsable=1470.0 

# températures en degrés
Tmoy=11
Tamp=9.1

# a diffusivité thermique du sous-sol
# a est exprimée en m2/s
a = lambda_sable/(cpsable*rho_sable)
# w pulsation exprimée en 1/s
w = 2*math.pi/(npy*step)

# temps correspondant au point le plus froid de l'année
tf = 18*24*step

titre="diffusivité : {} m2/s - pulsation : {} s-1".format("%.2E" % a,"%.2E" % w)

# profondeur équivalente liée à a ??
za = math.sqrt( 2 * a / w )


def Tsous_sol(z,i):
        
    """
    La température du sous-sol est une fonction sinusoidale :
    - de période = une année
    - de pulsation w ,obtenue à partir de l'équation de la propagation de la chaleur dans le sous-sol.
    
    z : profondeur en mètres
    
    i : indice temporel dans la discrétisation
        
    """
    
    factor = math.cos(w * (step * i - tf) + z/za ) * math.exp(z / za)
        
    return Tmoy - Tamp * factor 
        
    
pertes=np.zeros(npy)
z=1
for i in range(pertes.shape[0]):
    pertes[i]=Tsous_sol(z,i)
    
plt.subplot(111)
plt.title("température du sous-sol à {} m de profondeur\n{}".format(z,titre))
plt.xlabel("temps - une unité = {} s".format(step))
plt.plot(pertes)
plt.show()