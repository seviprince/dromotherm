import numpy as np

   
def lambda_eq(lambda1,lambda2,phi):
    
    """
    Cette fonction calcule la conductivité thermique équivalente d'un mélange de 2 matériaux
       
    """
    """
    Parameters
    ----------
    lambda1 : float
       Conductivité thermique du matériau 1 (W/m.K)
    lambda2 : float
       Conductivité thermique du matériau 2 (W/m.K)
       
    phi: float
        Porosité du mélange (-)
    """
    
     

    C=1.4

    B=C*(1-phi)**(10/9)/phi

    x=1-B*lambda1/lambda2

    y1=(1-lambda1/lambda2)*B*np.log(lambda2/(B*lambda1))/x**2

    y2=-(B+1)/2

    y3=-(B-1)/x

    coeff=(1/x)*2*(1-phi)**0.5

    return lambda1*(1-(1-phi)**0.5+coeff*(y1+y2+y3))

