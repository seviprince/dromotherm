# -*- coding: utf-8 -*-

# CODE DE RESOLUTION DU SYSTEME D'EQUATION DIFFERENTIELLE AU NIVEAU DE L'ECHANGEUR DROMOTHERM

from scipy.integrate import odeint 
import numpy as np
import matplotlib.pyplot as plt
import xlrd

# fonction définissant le système d'équations différentielles 
def F(X,temps):
    x,y,z = X[0],X[1],X[2] 
    
    return [(1/A1)*(rd_s*(y-x)-Hv*(x+273.15)+ks*(x-y)-epsilon*sigma*(x+273.15)**4+B1),-(1/B2)*(rd_s*(y-x)+rd_b*(y-z)-k*(y+273.15)/((rho_Cp_f*qf)/S+k/2)),-rd_b*(z-y)/C3]

# LES DONNEES NUMERIQUES
# LES DONNES CONSTANTES
rho_Cp_s=22546665 # capacité calorifique volumique de la couche de surface en (J/m^3.K)
rho_Cp_d=1769723 # capacité calorifique volumique de la couche drainante en (J/m^3.K)
rho_Cp_b=2676723 # capacité calorifique volumique de la couche de base en (J/m^3.K)
rho_Cp_f=4181000 # capacité calorifique volumique de couche deen (J/m^3.K)
h_s=0.06 # épaisseur de la couche de surface en m
h_d= 0.08 # épaisseur de la couche drainante en m
h_b=1 # épaisseur de la couche de base
albedo=0.08 # pas d'unité
epsilon=0.92 # émissivité (pas d'unité)
sigma=5.67*10**(-8) # constante de Stephan-Boltzmann
ks=2.34 #  conductivité thermique de la couche de surface en W/m.K
kd=1.56 #  conductivité thermique de la couche de drainante en W/m.K  
kb=1.76 #  conductivité thermique de la couche de base en W/m.K
qf=50/3600 # en m^3/s
k=10.83 # coefficient d'échange entre le fluide et la couche drainante 
mf=qf # débit massique mf =1kg/l * qf 
rd_s=27 # coefficient d'échange surfacique entre la couche drainate et la surface
rd_b=3.23 # coefficient d'échange surfacique entre la couche drainate et la surface
S=4 # surface de chaussée

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
B1=[] # B1 est créé pour des calculs intermédiaires pour des raisons de clarté
Hv=[] # Coefficients convectifs entre la surface et l'air: il varie en fction de la vitesse du vent
B1=np.asarray(B1)
Hv=np.asarray(Hv)
Hv=5.8+4.1*Vent
B1=Hv*(theta_air+273.15)+Rat+(1-albedo)*Rg 

# Calcul des solutions
solution = odeint(F,[5,5,5], temps)
temps=temps/3600 # on ramène le temps en heure

#Graphique
plt.xlabel('Temps (en seconde)')
plt.ylabel('Températures (en °C)')
plt.title("Profils de températures")
plt.plot(temps,solution)
plt.show()




