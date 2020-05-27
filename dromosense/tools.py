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
    
    etc etc
    
    return : vecteur numpy du besoin instantanné de chauffage en W
    
    """
    Rthe=1/(1/(Rm+Ri)+1/Rf)# Résistance thermique équivalente 
    return (Tconsigne-Text)/Rthe

def sol_tridiag(A,B,C,D):
    """
    Résout un système matriciel tridiagonal

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