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
