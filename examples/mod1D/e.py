from dromosense.tools import *
dt=3600
nt=1000
L=4
dx=0.75
qf=0.035/3600
dromo=OneDModel('meteo.txt',dt,nt,L,dx,qf)
Y=dromo.f1
dromo.f2
K=dromo.T[0,:,:]

K1=dromo.iterate(900,10)

