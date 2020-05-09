from dromosense import getCsvDatas

step, datas = getCsvDatas("meteo.csv",preview=True)

"""
col 0 : Température air (°C)
col 1 : Température point de rosée (°C)
col 2 : Nature des précipitations
col 3 : Vitesse du vent (m/s)
col 4 : Rayonnement global (W/m2)
col 5 : Rayonnement atmosphérique (W/m2)
"""

"""
Hv Coefficient convectif entre la surface et l'air (fonction affine de la vitesse du vent)
"""
Hv=5.8+4.1*datas[:,3]
B1=Hv*datas[:,0]+datas[:,5]+(1-albedo)*datas[:,4]

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
C_s=2144309
C_d=1769723
C_b=2676728
C_f=4181000

"""
épaisseurs en m
"""
h_s=0.06
h_d= 0.08
h_b=10

"""
albedo et émissivité
grandeurs sans unité
"""
albedo=0.08
epsilon=0.92

# constante de Stefan-Boltzmann en W/(m2K4)
#sigma=5.67*10**(-8)
sigma=5.67e-8

# température d'injection du fluide dans le dromotherme en °C
T_injection=10

# surface du dromotherme en m2
S_dromo=8
"""
surface d'échange entre le fluide et une paroi solide (couche drainante)
On suppose qu'elle est 100 fois plus grande de la surface du dromotherme et on arrondi !
"""
S=1000
# coefficient d'échange convectif entre le fluide et la couche drainante en W/(m2K)
h=2000