import numpy as np
import matplotlib.pyplot as plt
import xlrd

#IMPORTATION DES DONNES METEOS (VARIABLES EN FONCTION DU TEMPS)
fichier_excel = xlrd.open_workbook("meteo.xlsx") # appel du fichier meteo
feuille_1 = fichier_excel.sheet_by_index(0)
col=feuille_1.ncols # verification du nombre de colonnes et de lignes
ligne=feuille_1.nrows

# recuperation des differntes valeurs contenues dans chaque colonne
temps=[] # liste contenant les temps de mesures
theta_air=[] # liste contenant les temperatures de l'air
Vent=[] # vitesse du vent 
Rg=[] # Rayonnement solaire global
Rat=[] # Rayonnnement atmospherique 

for r in range(1, ligne):
    temps += [feuille_1.cell_value(rowx=r, colx=0)]
    theta_air+= [feuille_1.cell_value(rowx=r, colx=1)]
    Vent += [feuille_1.cell_value(rowx=r, colx=4)]
    Rg += [feuille_1.cell_value(rowx=r, colx=5)]
    Rat += [feuille_1.cell_value(rowx=r, colx=6)]
 
# conversion du type list en type array
temps=np.asarray(temps)
theta_air=np.asarray(theta_air)
Vent=np.asarray(Vent)
Rg=np.asarray(Rg)
Rat=np.asarray(Rat)

#AUTRES DONNEES NUMERIQUES

Hv=np.array([]) # Coefficients convectifs entre la surface et l'air: il varie en fction de la vitesse du vent
#Hv=5.8+4.1*Vent
Hv=4.65+3.96*Vent

rho_Cp_s=2144309.0 # capacite calorifique volumique de la couche de surface en (J/m^3.K)
rho_Cp_d=1769723.0 # capacite calorifique volumique de la couche drainante en (J/m^3.K)
rho_Cp_b=2676728.0 # capacite calorifique volumique de la couche de base en (J/m^3.K)
rho_Cp_f=4180000.0 # capacite calorifique volumique de couche du fluide (J/m^3.K)

L=4.0 # Largeur de chaussee
dx=0.4 # pas d'espace
h_s=0.06 # epaisseur de la couche de surface en m
h_d=0.08 # epaisseur de la couche drainante en m
h_b=10 # epaisseur de la couche de base

albedo=0.08 # pas d'unite
epsilon=0.92 # emissivite (pas d'unite)
sigma=5.67*10**(-8) # constante de Stefan-Boltzmann


ks=2.34 #  conductivite thermique de la couche de surface en W/m.K
kd=1.56 #  conductivite thermique de la couche de drainante en W/m.K  
kb=1.76 #  conductivite thermique de la couche de base en W/m.K

rd_s=2*ks*kd/(h_s*kd+h_d*ks) # coefficient d'echange surfacique entre la couche drainante et la surface

rd_b=2*kb*kd/(h_b*kd+h_d*kb) # coefficient d'echange surfacique entre la couche drainante et la surface 

rho_Cp_d=1.*1769723.0
rd_b=1
qf= 0.0/3600         # debit volumique du fluide (m^3/s)



phi=0.27     #( porosite de la couche drainante)
p=0.03       # ( la pente )
Kdr=qf/(h_d*p*L) # la permeabilite de la couche drainante


# CALCUL DU SCHEMA NUMERIQUE
dt=3600.0
V=Kdr*p*rho_Cp_f*h_d/(phi*rho_Cp_f*h_d+(1-phi)*rho_Cp_d*h_d) # la vitesse V figurant dans l'equation de transport 2
V2= V*dt # on calcule V*dt pour fixer le pas de temps tel que V*dt<=dx 

Nx = int(L/dx) # Nombre de blocs
Nt = int((temps[2926]-temps[0])/dt) # le nomnre d'iterations en temps
N=3*(Nx) + 3 # taille des matrices A,X,D,alpha et beta

Ts = np.zeros((Nt+1,Nx+1))
Tf = np.zeros((Nt+1,Nx+1))
Tb = np.zeros((Nt+1,Nx+1))
        
# Calcul des differents elements des vecteurs diagonaux
A=np.zeros((N))
B=np.zeros((N))
C=np.zeros((N))

alpha=np.zeros((N))
beta=np.zeros((N))
X=np.zeros((N))


a2=(1.0/dt)*((1-phi)*rho_Cp_d*h_d+phi*rho_Cp_f*h_d)+rd_s+rd_b
a3=rho_Cp_b*h_b/dt+rd_b
B[0] = -rd_s
A[1] = 10000.0
B[1] = 0.0
C[1] = 0.0
A[2]=a3
C[2]=-rd_b

for j in range (3,N):
    if(j%3==0):
        B[j]=-rd_s
    if(j%3==1):
        A[j]=a2
        B[j]=-rd_b
        C[j]=-rd_s
    if(j%3==2):
        A[j]=a3
        C[j]=-rd_b
       
# Conditions initiales

for j in range(0,Nx+1):
    Ts[0][j] = 10 + 273.15
    Tf[0][j] = 10 + 273.15
    Tb[0][j] = 10 + 273.15

# Conditions aux limites
for n in range(0,Nt+1):
    Ts[n][0] = 10  + 273.15
    Tf[n][0] = 10  + 273.15
    Tb[n][0] = 10  + 273.15


for n in range(0,Nt):
    D=np.array([])
    j = 0
    A[3*j]=rho_Cp_s*h_s/dt+rd_s+Hv[n+1]+4*epsilon*sigma*(Ts[n][j])**3
    Ds=(1.0-albedo)*Rg[n+1]+Rat[n+1]+Hv[n+1]*(theta_air[n+1]+273.15)+3*epsilon*sigma*(Ts[n][j])**4 + rho_Cp_s*h_s*Ts[n][j]/dt
    Df= 10000.0*Tf[n+1][j]
    Db=(rho_Cp_b*h_b/dt)*(Tb[n][j])
    D=np.append(D,[Ds,Df,Db])
    for j in range (1,Nx+1):
        A[3*j]=rho_Cp_s*h_s/dt+rd_s+Hv[n+1]+4*epsilon*sigma*(Ts[n][j])**3
        Ds=(1-albedo)*Rg[n+1]+Rat[n+1]+Hv[n+1]*(theta_air[n+1]+273.15)+3*epsilon*sigma*(Ts[n][j])**4 + rho_Cp_s*h_s*Ts[n][j]/dt
        #Df=((1/dt)*((1-phi)*rho_Cp_d*h_d+phi*rho_Cp_f*h_d))*(Tf[n][j])+(Kdr*p*rho_Cp_f*h_d/dx)*(Tf[n][j-1])
        #Df=((1.0/dt)*((1.0-phi)*rho_Cp_d*h_d+phi*rho_Cp_f*h_d)-qf*rho_Cp_f*h_d/dx)*(Tf[n][j])+(qf*rho_Cp_f*h_d/dx)*(Tf[n][j-1])
        Df=((1.0/dt)*((1.0-phi)*rho_Cp_d*h_d+phi*rho_Cp_f*h_d)-qf*rho_Cp_f*1/(L*dx))*(Tf[n][j])+(qf*rho_Cp_f*1/(L*dx))*(Tf[n][j-1])
        #print Df
        Db=(rho_Cp_b*h_b/dt)*(Tb[n][j])
        D=np.append(D,[Ds,Df,Db])
        
    alpha[0]=A[0]
    beta[0]=D[0]/alpha[0]
    for i in range(0,N-1):
        alpha[i+1]=A[i+1]-(C[i+1]*B[i])/alpha[i]
        beta[i+1]=(D[i+1]-B[i]*beta[i])/alpha[i+1]
    
    X[N-1]=beta[N-1]
    for i in range(N-2,-1,-1):
        X[i]=beta[i]-(B[i]*X[i+1]/alpha[i])

    j = 0
    Ts[n+1][j]=X[j*3]
    Tf[n+1][j]=X[j*3+1]
    Tb[n+1][j]=X[j*3+2]

    for j in range (1,Nx+1):
        Ts[n+1][j]=X[j*3]
        Tf[n+1][j]=X[j*3+1]
        Tb[n+1][j]=X[j*3+2]


np.savetxt('Ts.txt', Ts-273.15, fmt='%.2e')
np.savetxt('Tf.txt', Tf-273.15, fmt='%.2e')
np.savetxt('Tb.txt', Tb-273.15, fmt='%.2e')
        
temps=temps/3600
 #Graphique
figure1 = plt.figure(figsize = (10, 5))
plt.xlabel('Temps (en heure)')
plt.ylabel('Températures (en °C)')
plt.title("Profils de températures ")

plt.plot(temps,Ts[:,Nx]-273.15,label="Température couche de surface") 
plt.plot(temps,Tf[:,Nx]-273.15,label="Température du fluide")
plt.plot(temps,Tb[:,Nx]-273.15,label="Température couche de base")
#plt.plot(temps,Q,label="Chaleur récupérée")
#plt.show()     
plt.legend(loc=1)  
   
P=rho_Cp_f*qf*(Tf[:,Nx]-Tf[:,0]) # Puissance instantanée     
Energie=(np.sum(P))*dt/(3600*1000)
print('Energie récupérée=', Energie,'kWh')

 #Graphique
#figure2 = plt.figure(figsize = (10, 5))
plt.xlabel('Temps (en heure)')
plt.ylabel('Puissance (en W)')
#♦plt.plot(temps,P,label="La puissanceTf")     
        





