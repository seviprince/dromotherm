import math
# diamètre en mm
D = 25E-3
# viscosité en m2/s
V = 1E-6

print("Nous avons des tubes de diamètre {} m".format(D))

def reynolds(debit):
    """
    Re=Vitesse*Diamètre/viscosité cinématique
    Vitesse=Débit/section
    Re=4*Débit/(pi*Diamètre*viscosité)
    """
    Re=4*debit/ (math.pi * D * V)
    print("pour un débit de {} m3/s, reynolds vaut {}".format(debit,Re))
    return Re

reynolds(7.5*0.035/3600)
reynolds(1.2*7.5*0.035/3600)
reynolds(1.2*4*0.035/3600)

