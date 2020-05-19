from dromosense import rd
from dromosense.constantes import *
import numpy as np
import matplotlib.pyplot as plt

_input = np.loadtxt('input.txt')
nc = _input.shape[0] # nombre de couches
ha = _input[:,0] # hauteur des couches
le = _input[:,1] # coef d'echanges des couches (derniere valeur non utilisee)
rc = _input[:,2] # capacites calorifiques des couches
print(ha)
print(le)
print(rd(ks,kd,ha[0],ha[1]))
print(rd(kd,kb,ha[1],ha[2]))
print(rd(kb,kb,ha[2],ha[3]))
