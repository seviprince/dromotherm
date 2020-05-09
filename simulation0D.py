# -*- coding: utf-8 -*-

# CODE DE RESOLUTION DU SYSTEME D'EQUATION DIFFERENTIELLE AU NIVEAU DE L'ECHANGEUR DROMOTHERM

from scipy.integrate import odeint 
import numpy as np
import matplotlib.pyplot as plt
import xlrd

# fonction définissant le système d'équations différentielles 
def F(X,temps):
    Tentree_fluide=10
    Surface=8
    x,y,z = X[0],X[1],X[2] 
    y1=(1/A1)*(rd_s*(y-x)-Hv[int(temps/3600)]*(x+273.15)-epsilon*sigma*(x+273.15)**4+B1[int(temps/3600)])
    #y2=-(1/B2)*(rd_s*(y-x)+rd_b*(y-z)-k*(y+273.15)/((rho_Cp_f*qf)/S+k/2))
    y2=-(1/(B2*Surface)*(rd_s*Surface*(y-x)+rd_b*Surface*(y-z)-(qf*rho_Cp_f*(Tentree_fluide-(y-(y-Tentree_fluide)*np.exp(-k*S/(qf*rho_Cp_f)))))))
    y3=-rd_b*(z-y)/C3
    return [y1,y2,y3]

# LES DONNEES NUMERIQUES
# LES DONNES CONSTANTES
    
rho_Cp_s=2144309 # capacité calorifique volumique de la couche de surface en (J/m^3.K)
rho_Cp_d=1769723 # capacité calorifique volumique de la couche drainante en (J/m^3.K)
rho_Cp_b=2676728 # capacité calorifique volumique de la couche de base en (J/m^3.K)
rho_Cp_f=4181000 # capacité calorifique volumique de couche du fluide (J/m^3.K)

h_s=0.06 # épaisseur de la couche de surface en m
h_d= 0.08 # épaisseur de la couche drainante en m
h_b=10 # épaisseur de la couche de base

albedo=0.08 # pas d'unité
epsilon=0.92 # émissivité (pas d'unité)
sigma=5.67*10**(-8) # constante de Stefan-Boltzmann

ks=2.34 #  conductivité thermique de la couche de surface en W/m.K
kd=1.56 #  conductivité thermique de la couche de drainante en W/m.K  
kb=1.76 #  conductivité thermique de la couche de base en W/m.K
k=2000 # coefficient d'échange entre le fluide et la couche drainante

qf=0.05/3600 # en m^3/s
mf=1*qf # débit massique mf =1kg/l * qf 

rd_s=27 # coefficient d'échange surfacique entre la couche drainante et la surface
rd_b=3.28 # coefficient d'échange surfacique entre la couche drainante et la surface
S=1000 # surface d'échange fluide/solide

A1=rho_Cp_s * h_s # A1,B2, C3 sont des valeurs intermédaires de calculs
B2=rho_Cp_d * h_d
C3=rho_Cp_b * h_b

#IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
fichier_excel = xlrd.open_workbook("meteo.xlsx") # appel du fichier météo
feuille_1 = fichier_excel.sheet_by_index(0)
col=feuille_1.ncols # vérification du nombre de colonnes et de lignes
ligne=feuille_1.nrows

# récupération des différntes valeurs contenues dans chaque colonne
temps=[] # liste contenant les temps de mesures
theta_air=[] # liste contenant les températures de l'air
Vent=[] # vitesse vent 
Rg=[] # Rayonnement solaire global
Rat=[] # Rayonnnement atmosphérique 

for r in range(1, ligne):
    temps += [feuille_1.cell_value(rowx=r, colx=0)]
    theta_air+= [feuille_1.cell_value(rowx=r, colx=1)]
    Vent += [feuille_1.cell_value(rowx=r, colx=4)]
    Rg += [feuille_1.cell_value(rowx=r, colx=5)]
    Rat += [feuille_1.cell_value(rowx=r, colx=6)]
 
# conversion du type list en type arryay
temps=np.asarray(temps)
theta_air=np.asarray(theta_air)
Vent=np.asarray(Vent)
Rg=np.asarray(Rg)
Rat=np.asarray(Rat)

# B1 est créé pour des calculs intermédiaires pour des raisons de clarté
# Coefficients convectifs entre la surface et l'air: il varie en fction de la vitesse du vent
Hv=5.8+4.1*Vent
B1=Hv*(theta_air+273.15)+Rat+(1-albedo)*Rg 

# Calcul des solutions
solution = odeint(F,[10,10,10], temps)
Tentree_fluide=10
theta_d=solution[:,1]
#Tsortie_fluide=Tentree_fluide+ k*(theta_d+273.15)/((rho_Cp_f*qf)/S+k/2)
Tsortie_fluide=theta_d-(theta_d-Tentree_fluide)*np.exp(-k*S/(qf*rho_Cp_f))
temps=temps/3600 # on ramène le temps en heure

#Graphique
figure = plt.figure(figsize = (10, 5))
plt.xlabel('Temps (en heure)')
plt.ylabel('Températures (en °C)')
plt.title("Profils de températures 0D")
plt.plot(temps,solution[:,0],label="Température couche de surface")
plt.plot(temps,solution[:,1],label="Température couche drainante")
plt.plot(temps,solution[:,2],label="Température couche de base")
plt.plot(temps,Tsortie_fluide,label="Température de sortie du fluide")
plt.legend(loc="upper right")
plt.show()





