import numpy as np
import matplotlib.pyplot as plt

from dromosense.constantes import *

from datetime import datetime
import time
from dateutil import tz
CET=tz.gettz('Europe/Paris')

"""
fonctions généralistes
"""
def tsToTuple(ts):
    """
    ts : unix time stamp en s

    return date tuple tm_year, tm_mon, tm_mday, tm_hour, tm_min, tm_sec, tm_wday, tm_yday, tm_isdst
    """
    _time=datetime.fromtimestamp(ts,CET)
    _tuple=_time.timetuple()
    return(_tuple)

def basicAgenda(nbpts,step,start,summerStart,summerEnd,schedule=np.array([[8,17],[8,17],[8,17],[8,17],[8,17],[-1,-1],[-1,-1]])):
    """
    building an agenda indicating wether people are present or not

    nbpts : number of points the agenda will store

    step : time step in seconds

    start : starting unix time stamp in seconds

    summerStart,summerEnd : unix time stamps to define the summer period

    schedule : numpy array of size (7,2) with the presence hours for each day of the week

    returns : numpy vector of size nbpoints with activity indication (1: human presence , 0: no presence)
    
    utilisation :
    ```
    start = 1483232400
    
    summerStart = 1496278800
    
    summerEnd=1504141200
    
    step=3600
    
    nbpts=365*24 # un an avec une donnée à pas horaire
    
    schedule=np.array([[6,17],[8,18],[8,17],[8,17],[8,17],[-1,-1],[-1,-1]])
    
    agenda=basicAgenda(nbpts,step,start,summerStart,summerEnd,schedule=schedule)
    ```

    """
    verbose=False

    agenda=np.zeros(nbpts)
    time=start
    tpl=tsToTuple(time)
    work=0

    # do we have a summer ?
    summer=False
    if start<summerStart<summerEnd<=start+nbpts*step:
        summer=True
    print(summer)

    # fetching which day of the week have no presence at all if any
    weekend=[]
    for i in range(schedule.shape[0]):
        if -1 in schedule[i]:
            weekend.append(i)

    if verbose:
        print(weekend)

    # initial condition
    horaires=schedule[tpl.tm_wday]
    if tpl.tm_hour in range(horaires[0],horaires[1]):
        if tpl.tm_wday not in weekend:
            work=1

    agenda[0]=work

    previous=tpl

    for i in range(1,nbpts):

        goToNexti=False

        if summer and time<=summerEnd and time>=summerStart:
            agenda[i]=0
            goToNexti=True

        if not goToNexti:
            tpl=tsToTuple(time)
            horaires=schedule[tpl.tm_wday]
            if verbose:
                print("we are day {}".format(tpl.tm_wday))
                print("{} vs {} and {} vs {}".format(tpl.tm_hour,horaires[1],previous.tm_hour,horaires[1]-1))
            if tpl.tm_hour==horaires[1] and previous.tm_hour==horaires[1]-1:
                if tpl.tm_wday not in weekend:
                    work=0
            if tpl.tm_hour==horaires[0] and previous.tm_hour==horaires[0]-1:
                if tpl.tm_wday not in weekend:
                    work=1
            agenda[i]=work
            previous=tpl

        if verbose:
            print(agenda[i])
            input("press a key")

        time+=step

    return agenda

"""
fonctions techniques
"""
def rd(k1,k2,h1,h2):
    """
    calcule le coefficient d'échange surfacique entre 2 couches de conductivité k1 et k2 et d'épaisseurs h1 et h2

    W/(m2K)
    """
    return 2*k1*k2/(h1*k2+h2*k1)

def besoin_bat(Tconsigne,Text,Rm,Ri,Rf):
    """
    Calcule les besoins du bâtiment avec le modèle RC

    Tconsigne : température de consigne en °C

    Text : vecteur numpy de la température extérieure

    Rm : Résistance thermique des murs (K/W)

    Ri : Résistance superficielle intérieure (K/W)

    Rf : résistance de fuite (infiltrations+vitre+renouvellement d'air) K/W

    return : vecteur numpy du besoin instantanné de chauffage en W

    Par analogie électrique, on assimile les températures à des tensions et les puissances à des intensités

    en première approximation, on a donc (Tint-Text)/(Rm+Ri) + (Tc-Text)/Rf + C dTint/dt = Qchauffage

    soit C dTint/dt = Qchauffage - (Tint-Text) * (1/(Rm+Ri) + 1/Rf)

    Pour maintenir Tint constante et égale à Tconsigne, on doit donc développer :

    Qchauffage = (Tconsigne-Text) * (1/(Rm+Ri) + 1/Rf)

    """
    return (Tconsigne-Text)*(1/(Rm+Ri)+1/Rf)

def sol_tridiag(A,B,C,D):
    """
    Résout un système matriciel de la forme MX=D avec M une matrice tridiagonale ayant:

    A: vecteur constituant la diagonale principale

    B: vecteur constituant la diagonale supérieure

    C: le vecteur constituant la diagonale inférieure

    """
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

def Tsorties_echangeur(Te1,Te2,mf1,mf2,Cp1,Cp2,eff):
    """
    Calcul les températures au niveau des sorties d'un échangeur thermique'

    Parameters
    ----------
    Te1 : Température d'entrée du fluide chaud

    Te2 : Température d'entrée du fluide froid

    mf1 : Débit massique du fluide chaud

    mf2 : Débit massique du fluide froid

    Cp1 : Capacité calorifique massique du fluide chaud

    Cp2 : Capacité calorifique massique du fluide froid

    eff : efficacité de l'échangeur


    Returns
    -------
    Ts1 : Température de sortie du fluide chaud

    Ts2 : Température de sortie du fluide froid


    """
    if (mf1*Cp1)<=(mf2*Cp2):
        Ts1=Te1-eff*(Te1-Te2)
        Ts2=Te2+(mf1*Cp1/(mf2*Cp2))*(Te1-Ts1)
    else:
        Ts2=Te2+eff*(Te1-Te2)
        Ts1=Te1+(mf2*Cp2/(mf1*Cp1))*(Te2-Ts2)
    return Ts1,Ts2

class OneDModel:
    """
    dromotherm 1D model
    """
    def __init__(self,fname,dt,nt,L=4,l=5,dx=0.75,qf=0.035*5/3600):
        """
        fname : nom du fichier contenant les paramètres définissant la chaussée.
        chaque ligne est une couche et contient 3 valeurs séparées par des espaces :
        hauteur, coeff d'échanges, capacité calo

        dt : pas de temps en secondes

        nt : nombre de points dans la discrétisation temporelle

        L : largeur de chaussée en m
        
        l: la longueur de la chaussée en m ; on prend l=5m pour avoir une surface de 20m2 capable de chauffer le bâtiment

        dx : pas d'espace en m

        qf : débit volumique du fluide (m^3/s)

        objets construits lors de l'initialisation :

        ha : vecteur des hauteurs des couches en m

        le : vecteur des coef d'echanges des couches (derniere valeur non utilisee) en W/(m2.K)

        rc : vecteur des capacités calorifiques des couches en J/(m3.K)

        A,B,C : vecteurs définissant le système tridiagonal

        A : diagonale

        B : sur-diagonale

        C : sous-diagonale

        T : tenseur du champ de température exprimé en Kelvin !!

        - axe 0 : Temps

        - axe 1 : nombre de couches ou z

        - axe 2 : nombre de blocs ou x

        """
        _input = np.loadtxt(fname)
        nx = int(L/dx)
        nc = _input.shape[0]

        self.f1 = np.zeros((nt))
        self.f2 = np.zeros((nt))

        self.T = np.zeros((nt,nc,nx)) + kelvin
        self.dt = dt
        self.L = L
        self.dx = dx
        self.qf = qf
        self.l=l
        
        self.ha = _input[:,0]
        self.le = _input[:,1]
        self.rc = _input[:,2]

        self.A = np.zeros((nc))
        self.B = np.zeros((nc))
        self.C = np.zeros((nc))

        self.A[1:nc-1] = dt * (self.le[0:nc-2] + self.le[1:nc-1])
        self.A[nc-1] = dt * self.le[nc-2]
        self.A = self.A + self.ha*self.rc
        self.B[0:nc-1] = - dt * self.le[0:nc-1]
        self.C[1:nc] = - dt * self.le[0:nc-1]

    def iterate(self,n,Tinj):
        """
        n : time index (number of time steps)

        Tinj : injection temperature expressed in K
        """
        lambd = 10000000 # raideur pour imposer Tinj
        nx = self.T.shape[2]
        dt = self.dt
        for j in range(0,nx):
            self.A[0] = dt * (self.f2[n] + self.le[0] + 4.0*epsilon*sigma*self.T[n-1,0,j]**3) + self.ha[0] * self.rc[0]
            R = self.ha*self.rc*self.T[n-1,:,j]
            R[0] = R[0] + dt * (self.f1[n] + 3.0*epsilon*sigma*self.T[n-1,0,j]**4)
            if j==0:
               R[1] = dt*lambd*Tinj
               self.A[1] = dt*lambd*1.0
               self.C[1] = 0.0
               self.B[1] = 0.0
            else:
               R[1] = R[1] + dt * (self.qf * Cpf * rho_eau /(self.l*self.dx)) * (self.T[n-1,1,j-1]-self.T[n-1,1,j])
               self.C[1] = - dt * self.le[0]
               self.B[1] = - dt * self.le[1]
               self.A[1] = dt * (self.le[0] + self.le[1]) + self.ha[1] * self.rc[1]
            self.T[n,:,j] = sol_tridiag(self.A,self.B,self.C,R)

    def showVectors(self):
        """
        affiche les vecteurs A,B,C
        """
        print("A is {}".format(self.A))
        print("B is {}".format(self.B))
        print("C is {}".format(self.C))
