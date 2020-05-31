import numpy as np
import matplotlib.pyplot as plt

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
        Ts2=Te2-eff*(Te1-Te2)
        Ts1=Te1+(mf2*Cp2/(mf1*Cp1))*(Te1-Ts1)
    return Ts1,Ts2    

