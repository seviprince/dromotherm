import numpy as np
from dromosense.constantes import kelvin, Z0,Zbe_i,Zi_1,Z1_2,Z2_3,Z3_4,Z4_5,Z5_i,Zs_sol,Tf, e_c,L,Cpf
import math
from scipy.optimize import fsolve


class SenCityTwo:
    """
    Cette classe réalise un DEUXIEME couplage du système complet et nourrit un ensemble de vecteurs numpy:
    
    Tinj_sto : température du fluide de charge à l'injection dans le stockage
    
    Tsor_sto : température du fluide de charge après transit dans le stockage
    
    Tsor_pac : température du fluide de décharge au sortir de la PAC à la réinjection dans le stockage
    
    Tinj_pac : température du fluide de décharge à l'injection dans la PAC
    
 
    Tinj_dro : température d'injection dans le dromotherme
    
    Tsor_dro : température de sortie du dromotherme
    
    Ts: température de la partie englacée de chaque couche du stockage
    
    Tl: Température de la partie non englacée de chaque couche du stcoakge
    
    zfu: position du front de glace dans chaque couche du stockage
    
    run_z1: booléen qui nous dit s'il y a apparition ou pas du gel à l'interface des couches 1 et 2

    run_z3: booléen qui nous dit s'il y a apparition ou pas du gel à l'interface des couches 3 et 4
    
    On initialise Tinj_dro et Tsor_dro à 10
    
    pertes : pertes instantannées dans le stockage en W (prérequis : nécessite l'utilisation de setPertes)
    
    agenda_dro et agenda_pac : agendas de fonctionnement du dromotherme et de la PAC
    
    Pgeo : puissance géothermique à extraire pour satisfaire le besoin du bâtiment en W
    

    
    Pour l'utiliser :
    
    ```
    # instanciation
    # RSB = route stock bâtiment
    RSB=SenCityTwo(3600,meteo.shape[0],nc)
    # définition du système
    RSB.set(eff,k,coeff,cste,lambda_be,ath_be,lambda_iso,ath_iso,lambda_l,ath_l,rho_l,Cp_l,lambda_s,ath_s,rho_s,msto,cpsto,mpac,cpac,Rthf,Sb)
    # injection du besoin
    RSB.Pgeo = Pgeo (via les données expérimentales)
    # définition de la fenêtre de simulation et paramétrage des agendas
    simStart = i_summerStart
    simEnd=i_summerStart+365*24
    RSB.agenda_dro[simStart:simEnd]=np.ones(simEnd-simStart)
    RSB.agenda_pac[simStart:simEnd]=np.ones(simEnd-simStart)
    for i in range(simStart,simEnd):
        if RSB.Pgeo[i]==0:
            RSB.agenda_pac[i]=0    
    # bouclage
    # dromo : échangeur dromotherme selon la classe OneDModel, traversé par un débit unitaire qdro_u en m3/s
    # unitaire s'entendant par mètre linéaire selon le profil en long
    for i in range (int(simStart),int(simEnd)):
        RSB.SystemLoop(n,qdro_u,dromo)
    ```
    """    

            
    def __init__(self,step,size,nc):
        
        self.size=size
        self.step=step
        self.nc=nc
        self.Tinj_sto=5*np.ones((size,2))
        self.Tsor_sto=5*np.ones((size,2))
        self.Tmoy_inj_sto=5*np.ones(size)
        self.Tmoy_sor_sto=5*np.ones(size)
        
        self.Tsor_pac=-1*np.ones((size,2))
        self.Tinj_pac=-1*np.ones((size,2))
        self.Tmoy_inj_pac=-1*np.ones(size)
        self.Tmoy_sor_pac=-1*np.ones(size)       
        self.Tinj_dro=10*np.ones(size)
        self.Tsor_dro=10*np.ones(size)  
        
        self.Ts=0*np.ones((size,nc))
        self.Tl=5*np.ones((size,nc))
        self.Tl[:,0]=6
        self.T_inter=5*np.ones((size,7))
        self.T_be=5*np.ones(size)
        self.T_iso1=5*np.ones(size)
        self.T_iso2=5*np.ones(size)
        
        
        self.Text=5*np.ones(size)
        

        self.agenda_dro=np.zeros(size)
        self.agenda_pac=np.zeros(size)
        
        self.zfu=np.zeros((size,nc))
        self.zfu[:,0]=Z1_2
        self.zfu[:,1]=Z1_2
        self.zfu[:,2]=Z3_4
        self.zfu[:,3]=Z3_4
        self.zfu[:,4]=Z4_5

        self.run_z1=True
        self.run_z3=True
        
        self.Pgeo=0*np.ones(size)
        
        #self.Pgeo[int(self.size/3):int(self.size)]=0
        self.P_1=1000*np.ones(size)        
        self.P_2=1000*np.ones(size)



    def set(self,eff,k,coeff,cste,lambda_be,ath_be,lambda_iso,ath_iso,lambda_l,ath_l,rho_l,Cp_l\
        ,lambda_s,ath_s,rho_s,msto,cpsto,mpac,cpac,Rthf,Sb):
        
        """
        paramètres permettant de caractériser le système couplé
        
        eff : efficacité de l'échangeur de séparation de réseaux
        
        k : coefficient du système géothermique équipant le stockage en W/K
        
        coeff : sans unité - rapport des conductivités thermiques des fluides circulant dans l'échangeur de séparation de réseaux (dromo/stockage)
        
        msto : débit massique du fluide dans le système géothermique du stockage en kg/s
        
        cpsto : capacité thermique du fluide circulant dans le système géothermique du stockage en J/(K.kg)
        
        msable : masse de sable humide (sable sec + eau ) en kg
        
        cpsable : capacité calorifique massique du sable humide en J/(K/kg)
        
        mpac : débit massique du fluide dans la PAC en kg/s
        
        cpac : capacité calorifique massique du fluide dans la PAC en J/(K.kg)
        
        """
        
        self.eff=eff
        self.k=k
        self.coeff=coeff
        self.cste=cste
        self.lambda_be=lambda_be
        self.ath_be=ath_be
        self.lambda_iso=lambda_iso
        self.ath_iso=ath_iso
        self.lambda_l=lambda_l
        self.ath_l=ath_l
        self.rho_l=rho_l
        self.lambda_s=lambda_s
        self.ath_s=ath_s
        self.rho_s=rho_s
        self.msto=msto
        self.cpsto=cpsto
     
        self.mpac=mpac
        self.cpac=cpac
        
        self.Rthf=Rthf
        self.Sb=Sb
        
        self.Cp_l=Cp_l
        
        self.Cp_f=cpac
       
        
        """
        calcul de B et C - cf notebook
        """
        
        self.B = (msto * cpsto - k/2) * coeff * eff
        
        self.C = 1 - k / (2 * mpac *cpac)
        
        print("le k du système géothermique vaut {} W/K".format(k))
        print("coeff vaut {}".format(coeff))
        print("B vaut {} W/K".format(self.B))
        print("C vaut {}".format(self.C))        

        
  #----------------------------------------------------------------------------------------------------------------------  
    def coeffs(self,T,Tsup,Tinf,Zsup,Zinf):
        
        """
        La température de chaque couche est sous une forme polynomiale T=az^2+b*z+c
        Cette fonction  calcule les coefficients a, b et c en résolvant le système d'équations:
           
            | a*Zsup^2+b*Zsup+c=Tsup
            | a*Zinf^2+b*Zinf+c=Tinf
            | a*(Zsup^2+Zsup*Zinf+Zinf^2)*3+b*(Zsup+Zinf)/2=T

    
        Parameters
        ----------
        T : float
            Température de la couche (en K)
          
        Tsup : float
            Température de l'interface supérieure de la couche
         
        Tinf : float
            Température de l'interface inférieure de la couche 
         
        Zsup :float
            Position ou côte de l'interface supérieure de la couche
           
        Zinf :float
            Position ou côte de l'interface inférieure de la couche
           
        """
        
        self.T=T
        self.Tsup=Tsup
        self.Tinf=Tinf
        self.Zsup=Zsup
        self.Zinf=Zinf
        
        

    
        if self.Zinf==0.0 or self.Zsup==0.0:
            self.a=3*(-2*self.T+self.Tsup+self.Tinf)/(self.Zsup-self.Zinf)**2
            self.b=2*(3*self.T-self.Tsup-2*self.Tinf)/(self.Zsup-self.Zinf)
            self.c=self.Tsup-self.a*self.Zsup**2-self.b*self.Zsup
   
        else:
    
   
            x1=Zinf**2
            x2=Zsup**2
            x3=(Zinf**2+Zinf*Zsup+Zsup**2)/3
            
            y1=Zinf
            y2=Zsup
            y3=(Zinf+Zsup)/2
        
            Y2=y2*x1-y1*x2
            W2=x1-x2
            Y3=x1*y3-x3*y1
            W3=x1-x3
            
            Z2=x1*Tsup-x2*Tinf
            Z3=x1*T-x3*Tinf
    
            self.c=(Y2*Z3-Y3*Z2)/(Y2*W3-Y3*W2-1e-30) # j'ai fait -1-e-15 pour eviter les warnings du type : divide by zero
            self.b=(Z2-W2*self.c)/Y2
            self.a=(Tinf-y1*self.b-self.c)/x1
            # On peut aussi utiliser la fonction linalg de numpy mais elle est coûteuse en temps d'éxécution     
        return self.a,self.b,self.c
  #----------------------------------------------------------------------------------------------------------------------  
    
    def T_interface(self,Tsup1,T1,Tinter,T2,Tinf2,Zsup1,Zinter,Zinf2,lambda1,lambda2,**kwargs):
   
        """
        Cette fonction calcule la température du tube à l'interface de chaque couche
    
        Parameters
        ----------
        Tsup1 : float
            Température de l'interface supérieure de la couche 1 (°C)
            
        T1 : float
            Température de la couche 1 (°C)
            
        Tinter : float
            Température à l'interface de la couche 1 et de la couche 2 (°C)
        T2 : float
            Température de la couche 2  (°C)     
        Tinf2 : 
            Température de l'interface supérieure de la couche 1 (°C)
        Zsup1 : float
            Position de l'interface supérieure de la couche 1 (m)
        Zinter : float
            Position de l'interface des couches 1 et 2 (m)
        Zinf2 : float
            Position de l'interface inférieure de la couche 2(m)
        lambda1 : float
            Conductivité thermique de la couche 1 (W/m.K)
        lambda2 : float
            Conductivité thermique de la couche 2 (W/m.K)
        **kwargs : Paramètres facultatifs
            m: le débit massique du fluide (kg/s) 
            Cp: la capacité thermique massique du fluide (J/Kg.K)
   
        """
    

        self.Tsup1=Tsup1
        self.T1=T1
        self.Tinter=Tinter
        self.T2=T2
        self.Tinf2=Tinf2
        self.Zsup1=Zsup1
        self.Zinter=Zinter
        self.Zinf2=Zinf2
        self.lambda1=lambda1
        self.lambda2=lambda2
        
    
        self.a1,self.b1,self.c1=self.coeffs(T1,Tsup1,Tinter,Zsup1,Zinter)
        
        self.a2,self.b2,self.c2=self.coeffs(T2,Tinter,Tinf2,Zinter,Zinf2)
        
            
                   
        if "m" not in kwargs or "Tin" not in kwargs :
                           
            return self.lambda1*(2*self.a1*self.Zinter+self.b1)-self.lambda2*(2*self.a2*self.Zinter+self.b2)
    
        
        else:
            
            self.m=kwargs["m"]
            self.Tin=kwargs["Tin"]
    
            if self.m==0.0:
                return self.lambda1*(2*self.a1*self.Zinter+self.b1)-self.lambda2*(2*self.a2*self.Zinter+self.b2)
            else:
                
                A=(1/(2*self.m*Cpf)+self.Rthf)*self.Sb
                B=-self.lambda1*(2*self.a1*self.Zinter+self.b1)+self.lambda2*(2*self.a2*self.Zinter+self.b2)
               # print(B*self.Sb)
                return self.Tinter-self.Tin-A*B
#----------------------------------------------------------------------------------------------------------------------
    def TfluidePac(self,Tinter,Tin,m,Q):
        
       
        if m==0:
            Tout=Tinter
            Tin_new=Tinter
    
        else:
      
            A=m*self.cpac*self.Rthf-0.5
            B=m*self.cpac*self.Rthf+0.5
            Tout=(Tinter+A*Tin)/B                      
            Tin_new=Tout-Q/(2*m*self.cpac)
            
        return Tout, Tin_new
#----------------------------------------------------------------------------------------------------------------------    
    def TfluideSto(self,Tinter,Tin,Tsor_dro,m,eff):
    
       
        if m==0:
            Tout=Tinter
            Tin_new=Tinter
    
        else:
            
            A=m*self.Cp_f*self.Rthf-0.5
            B=m*self.Cp_f*self.Rthf+0.5     

            Tout=(Tinter+A*Tin)/B                      
            Tin_new=Tsor_dro-self.eff*(Tsor_dro-Tout)
            
        return Tout, Tin_new     
   
#----------------------------------------------------------------------------------------------------------------------
    def equations(self,Y,*p):
            
        
       """
        Définit le système des  équations décrivant les transferts thermiques dans le stockage avec ou
        
        sans changement de phase
    
        Y: vecteur contenant les 23 inconnues à l'instant n+1 (Tbe,Tbe_i,Tiso1,Ts1,Tl1,zfu1,....,Tiso2)
        
        Tout_pac=Tinj_pac
        Tin_pac=Tsor_pac
        
    
        """ 
        
        
       (Tbe,Tbe_i,Tiso1,Ti_1,Ts1,Tl1,zfu1,T1_2,Ts2,Tl2,zfu2,T2_3,Ts3,Tl3,zfu3,T3_4,Ts4,Tl4,zfu4,T4_5,Tl5\
       ,T5_i,Tiso2,Tinj_pac1,Tinj_pac2,Tsor_pac1,Tsor_pac2,Tinj_sto1,Tinj_sto2,Tsor_sto1,Tsor_sto2)=Y
       
       (n,p1,p2,p3,p4)=p
       
       Ts_sol=self.Tsous_sol(-Zs_sol,n+1)

       if p==(n,1,0,1,0):
           
           a_be,b_be,c_be=self.coeffs(Tbe,self.Text[n+1],Tbe_i,Z0,Zbe_i)
           a_iso1,b_iso1,c_iso1=self.coeffs(Tiso1,Tbe_i,Ti_1,Zbe_i,Zi_1)
                
           a_1,b_1,c_1=self.coeffs(Tl1,Ti_1,T1_2*p1+p2*Tf,Zi_1,Z1_2*p1+p2*zfu1)
           a_p1,b_p1,c_p1=(0,0,0)
             
           a_p2,b_p2,c_p2=(0,0,0)
           a_2,b_2,c_2=self.coeffs(Tl2,T1_2*p1+p2*Tf,T2_3,Z1_2*p1+p2*zfu2,Z2_3)
             
           a_3,b_3,c_3=self.coeffs(Tl3,T2_3,p3*T3_4+p4*Tf,Z2_3,p3*Z3_4+p4*zfu3)
           a_p3,b_p3,c_p3=(0,0,0)
             
           a_p4,b_p4,c_p4=(0,0,0)
           a_4,b_4,c_4=self.coeffs(Tl4,T3_4*p3+p4*Tf,T4_5,p3*Z3_4+p4*zfu4,Z4_5)
                
           a_5,b_5,c_5=self.coeffs(Tl5,T4_5,T5_i,Z4_5,Z5_i)
           a_iso2,b_iso2,c_iso2=self.coeffs(Tiso2,T5_i,self.Tsous_sol(-Zs_sol,n+1),Z5_i,Zs_sol) 
          # a_iso2,b_iso2,c_iso2=self.coeffs(Tiso2,T5_i,Ts_sol,Z5_i,Zs_sol) 
           

           
       if p==(n,0,1,1,0):
           a_be,b_be,c_be=self.coeffs(Tbe,self.Text[n+1],Tbe_i,Z0,Zbe_i)
           a_iso1,b_iso1,c_iso1=self.coeffs(Tiso1,Tbe_i,Ti_1,Zbe_i,Zi_1)
                
           a_1,b_1,c_1=self.coeffs(Tl1,Ti_1,T1_2*p1+p2*Tf,Zi_1,Z1_2*p1+p2*zfu1)
           a_p1,b_p1,c_p1=self.coeffs(Ts1,Tf,T1_2,zfu1,Z1_2)
             
           a_p2,b_p2,c_p2=self.coeffs(Ts2,T1_2,Tf,Z1_2,zfu2)
           a_2,b_2,c_2=self.coeffs(Tl2,T1_2*p1+p2*Tf,T2_3,Z1_2*p1+p2*zfu2,Z2_3)
             
           a_3,b_3,c_3=self.coeffs(Tl3,T2_3,p3*T3_4+p4*Tf,Z2_3,p3*Z3_4+p4*zfu3)
           a_p3,b_p3,c_p3=(0,0,0)
             
           a_p4,b_p4,c_p4=(0,0,0)
           a_4,b_4,c_4=self.coeffs(Tl4,T3_4*p3+p4*Tf,T4_5,p3*Z3_4+p4*zfu4,Z4_5)
                
           a_5,b_5,c_5=self.coeffs(Tl5,T4_5,T5_i,Z4_5,Z5_i)
           a_iso2,b_iso2,c_iso2=self.coeffs(Tiso2,T5_i,self.Tsous_sol(-Zs_sol,n+1),Z5_i,Zs_sol) 
       #    a_iso2,b_iso2,c_iso2=self.coeffs(Tiso2,T5_i,Ts_sol,Z5_i,Zs_sol) 
         
       if p==(n,1,0,0,1):
           
           a_be,b_be,c_be=self.coeffs(Tbe,self.Text[n+1],Tbe_i,Z0,Zbe_i)
           a_iso1,b_iso1,c_iso1=self.coeffs(Tiso1,Tbe_i,Ti_1,Zbe_i,Zi_1)
                
           a_1,b_1,c_1=self.coeffs(Tl1,Ti_1,T1_2*p1+p2*Tf,Zi_1,Z1_2*p1+p2*zfu1)
           a_p1,b_p1,c_p1=(0,0,0)
             
           a_p2,b_p2,c_p2=(0,0,0)
           a_2,b_2,c_2=self.coeffs(Tl2,T1_2*p1+p2*Tf,T2_3,Z1_2*p1+p2*zfu2,Z2_3)
             
           a_3,b_3,c_3=self.coeffs(Tl3,T2_3,p3*T3_4+p4*Tf,Z2_3,p3*Z3_4+p4*zfu3)
           a_p3,b_p3,c_p3=self.coeffs(Ts3,Tf,T3_4,zfu3,Z3_4)
             
           a_p4,b_p4,c_p4=self.coeffs(Ts4,T3_4,Tf,Z3_4,zfu4)
           a_4,b_4,c_4=self.coeffs(Tl4,T3_4*p3+p4*Tf,T4_5,p3*Z3_4+p4*zfu4,Z4_5)
                
           a_5,b_5,c_5=self.coeffs(Tl5,T4_5,T5_i,Z4_5,Z5_i)
           a_iso2,b_iso2,c_iso2=self.coeffs(Tiso2,T5_i,self.Tsous_sol(-Zs_sol,n+1),Z5_i,Zs_sol) 
          # a_iso2,b_iso2,c_iso2=self.coeffs(Tiso2,T5_i,Ts_sol,Z5_i,Zs_sol) 
                              
       if p==(n,0,1,0,1):            
           a_be,b_be,c_be=self.coeffs(Tbe,self.Text[n+1],Tbe_i,Z0,Zbe_i)
           a_iso1,b_iso1,c_iso1=self.coeffs(Tiso1,Tbe_i,Ti_1,Zbe_i,Zi_1)
                
           a_1,b_1,c_1=self.coeffs(Tl1,Ti_1,T1_2*p1+p2*Tf,Zi_1,Z1_2*p1+p2*zfu1)
           a_p1,b_p1,c_p1=self.coeffs(Ts1,Tf,T1_2,zfu1,Z1_2)
             
           a_p2,b_p2,c_p2=self.coeffs(Ts2,T1_2,Tf,Z1_2,zfu2)
           a_2,b_2,c_2=self.coeffs(Tl2,T1_2*p1+p2*Tf,T2_3,Z1_2*p1+p2*zfu2,Z2_3)
             
           a_3,b_3,c_3=self.coeffs(Tl3,T2_3,p3*T3_4+p4*Tf,Z2_3,p3*Z3_4+p4*zfu3)
           a_p3,b_p3,c_p3=self.coeffs(Ts3,Tf,T3_4,zfu3,Z3_4)
             
           a_p4,b_p4,c_p4=self.coeffs(Ts4,T3_4,Tf,Z3_4,zfu4)
           a_4,b_4,c_4=self.coeffs(Tl4,T3_4*p3+p4*Tf,T4_5,p3*Z3_4+p4*zfu4,Z4_5)
                
           a_5,b_5,c_5=self.coeffs(Tl5,T4_5,T5_i,Z4_5,Z5_i)
           a_iso2,b_iso2,c_iso2=self.coeffs(Tiso2,T5_i,self.Tsous_sol(-Zs_sol,n+1),Z5_i,Zs_sol) 
          # a_iso2,b_iso2,c_iso2=self.coeffs(Tiso2,T5_i,Ts_sol,Z5_i,Zs_sol) 
                   
       #équations
         
       self.eq_be=Tbe-self.T_be[n]-2*self.ath_be*a_be*self.step
       self.eq_bi=self.T_interface(self.Text[n+1],Tbe,Tbe_i,Tiso1,Ti_1,Z0,Zbe_i,Zi_1,self.lambda_be,self.lambda_iso)
       self.eq_iso1=Tiso1-self.T_iso1[n]-2*self.ath_iso*a_iso1*self.step
         
       self.eqi_1=self.T_interface(Tbe_i,Tiso1,Ti_1,Tl1,T1_2*p1+p2*Tf,Zbe_i,Zi_1,Z1_2*p1+p2*zfu1,self.lambda_iso,self.lambda_l)
         
       self.eq_l1=Tl1-self.Tl[n,0]-2*self.ath_l*a_1*self.step
       self.eq_s1=(Ts1-self.Ts[n,0])-2*self.ath_s*a_p1*self.step*p2
       self.eq_zfu1=(zfu1-self.zfu[n,0])-self.cste*(self.lambda_s*(2*a_p1*zfu1+b_p1)-self.lambda_l*(2*a_1*zfu1+b_1))*p2
         
       self.eq_12=self.T_interface(Ti_1*p1+p2*Tf,Tl1*p1+p2*Ts1, T1_2,p1*Tl2+p2*Ts2, p1*T2_3+p2*Tf, p1*Zi_1+p2*zfu1, Z1_2\
                , p1*Z2_3+p2*zfu2,p1*self.lambda_l+p2*self.lambda_s,p1*self.lambda_l+p2*self.lambda_s\
                ,m=self.mpac*self.agenda_pac[n+1],Tin=Tsor_pac1)
            
       self.eq_l2=Tl2-self.Tl[n,1]-2*self.ath_l*a_2*self.step
       self.eq_s2=(Ts2-self.Ts[n,1])-2*self.ath_s*a_p2*self.step*p2
       self.eq_zfu2=(zfu2-self.zfu[n,1])-self.cste*(self.lambda_s*(2*a_p2*zfu2+b_p2)-self.lambda_l*(2*a_2*zfu2+b_2))*p2  
           
       self.eq_23=self.T_interface(p1*T1_2+p2*Tf, Tl2, T2_3, Tl3, p3*T3_4+p4*Tf, p1*Z1_2+p2*zfu2, Z2_3, p3*Z3_4+p4*zfu3\
                , self.lambda_l, self.lambda_l,m=self.msto*self.agenda_dro[n+1],Tin=Tinj_sto1)
         
       self.eq_l3=Tl3-self.Tl[n,2]-2*self.ath_l*a_3*self.step
       self.eq_s3=(Ts3-self.Ts[n,2])-2*self.ath_s*a_p3*self.step*p4
       self.eq_zfu3=(zfu3-self.zfu[n,2])-self.cste*(self.lambda_s*(2*a_p3*zfu3+b_p3)-self.lambda_l*(2*a_3*zfu3+b_3))*p4
         
       self.eq_34=self.T_interface(p3*T2_3+p4*Tf, p3*Tl3+p4*Ts3, T3_4,p3*Tl4+p4*Ts4, p3*T4_5+p4*Tf, p3*Z2_3+p4*zfu3\
                , Z3_4, p3*Z4_5+p4*zfu4, p3*self.lambda_l+p4*self.lambda_s,  p3*self.lambda_l+p4*self.lambda_s\
                ,m=self.mpac*self.agenda_pac[n+1],Tin=Tsor_pac2)
         
       self.eq_l4=Tl4-self.Tl[n,3]-2*self.ath_l*a_4*self.step
       self.eq_s4=(Ts4-self.Ts[n,3])-2*self.ath_s*a_p4*self.step*p4
       self.eq_zfu4=(zfu4-self.zfu[n,3])-self.cste*(self.lambda_s*(2*a_p4*zfu4+b_p4)-self.lambda_l*(2*a_4*zfu4+b_4))*p4   
         
       self.eq_45=self.T_interface(p3*T3_4+p4*Tf, Tl4, T4_5, Tl5, T5_i, p3*Z3_4+p4*zfu4, Z4_5, Z5_i, self.lambda_l, self.lambda_l\
                ,m=self.msto*self.agenda_dro[n+1],Tin=Tinj_sto2)
         
       self.eq_l5=Tl5-self.Tl[n,4]-2*self.ath_l*a_5*self.step

         
       self.eq_5i=self.T_interface(T4_5,Tl5,T5_i,Tiso2,Ts_sol,Z4_5,Z5_i,Zs_sol,self.lambda_l,self.lambda_iso)
       self.eq_5i=self.T_interface(T4_5,Tl5,T5_i,Tiso2,Ts_sol,Z4_5,Z5_i,Zs_sol,self.lambda_l,self.lambda_iso)
         
       # isolant 2
         
       self.eq_iso2=Tiso2-self.T_iso2[n]-2*self.ath_iso*a_iso2*self.step
       
       # Les fluides
       

          
       Tpo1,Tpi1=self.TfluidePac(T1_2,Tsor_pac1,self.mpac*self.agenda_pac[n+1],self.Pgeo[n+1])
       Tpo2,Tpi2=self.TfluidePac(T3_4,Tsor_pac2,self.mpac*self.agenda_pac[n+1],self.Pgeo[n+1])
       
       self.eq_inj_pac1=Tinj_pac1-Tpo1
       self.eq_inj_pac2=Tinj_pac2-Tpo2
       self.eq_sor_pac1=Tsor_pac1-Tpi1
       self.eq_sor_pac2=Tsor_pac2-Tpi2
    
      
       Tdo1,Tdi1=self.TfluideSto(T2_3,Tinj_sto1,self.Tsor_dro[n+1],self.msto*self.agenda_dro[n+1],self.eff)
       Tdo2,Tdi2=self.TfluideSto(T4_5,Tinj_sto2,self.Tsor_dro[n+1],self.msto*self.agenda_dro[n+1],self.eff)

        
       self.eq_sor_sto1=Tsor_sto1-Tdo1
       self.eq_sor_sto2=Tsor_sto2-Tdo2     
       self.eq_inj_sto1=Tinj_sto1-Tdi1      
       self.eq_inj_sto2=Tinj_sto2-Tdi2

      # resultat=(self.eq_be,self.eq_bi,self.eq_iso1,self.eqi_1,self.eq_s1,self.eq_l1,self.eq_zfu1,self.eq_12,self.eq_s2,self.eq_l2,self.eq_zfu2,self.eq_23,self.eq_s3,self.eq_l3,self.eq_zfu3,self.eq_34,self.eq_s4,self.eq_l4,self.eq_zfu4,self.eq_45,self.eq_l5,self.eq_5i,self.eq_iso2, self.eq_inj_pac1, self.eq_inj_pac2, self.eq_sor_pac1, self.eq_sor_pac2,self.eq_inj_sto1, self.eq_inj_sto2,self.eq_sor_sto1,self.eq_sor_sto2)
      # res = [type(ele) for ele in resultat]
 
        # printing result
       #print("The data types of tuple in order are : " + str(res))
    

       return self.eq_be,self.eq_bi,self.eq_iso1,self.eqi_1,self.eq_s1,self.eq_l1,self.eq_zfu1,self.eq_12\
              ,self.eq_s2,self.eq_l2,self.eq_zfu2,self.eq_23,self.eq_s3,self.eq_l3,self.eq_zfu3,self.eq_34\
              ,self.eq_s4,self.eq_l4,self.eq_zfu4,self.eq_45,self.eq_l5,self.eq_5i,self.eq_iso2\
              ,self.eq_inj_pac1, self.eq_inj_pac2, self.eq_sor_pac1, self.eq_sor_pac2,self.eq_inj_sto1\
              ,self.eq_inj_sto2,self.eq_sor_sto1,self.eq_sor_sto2
    #----------------------------------------------------------------------------------------------------------------------------------
      


    def SystemLoop(self,n,qdro_u,dromo):
        """
        dromo : échangeur dromotherme simulé selon la classe OneDModel
        
        qdro_u : débit unitaire en m3/s traversant le dromotherme, unitaire s'entendant par mètre linéaire selon le profil en long
        
        1) On applique StockLoop avec les résutats de l'état précédant, ce qui nous permet de calculer Tsable[i]
        
        2) On met à jour les température d'injection et de sortie de la PAC
        
        3) On réalise ensuite une itération de dromotherme selon 2 cas distincts
        
        ********************************************************************************
        cas 1: le dromotherme est en marche
        
        le fluide circule avec un débit unitaire qdro_u
          
        Si Tsor_dro <= Tsable : goto cas 2 (dromotherme à l'arrêt) + on passe la valeur de agenda_dro[i] à 0
          
        Si Tsor_dro > Tsable : alimentation du stockage
        ********************************************************************************       
        cas 2: le dromotherme est à l'arrêt - débit nul
        
        l'échangeur de séparation de réseau ne tourne pas
        
        pas de prélèvement par l'échangeur de séparation de réseau
        
        `Tinj_dro[i] = Tsor_dro[i]`
        
        pas d'évolution des températures du stockage
        
        `Tsor_sto[i]=Tsor_sto[i-1]` et `Tinj_sto[i]=Tinj_sto[i-1]`
          
        """
        # étape 1 
     
        dro=self.agenda_dro[n]
        pac=self.agenda_pac[n]
 
        if self.T_inter[n,2]>=0 and self.T_inter[n,4]>=0:
            
            p1=1
            p2=0
            p3=1
            p4=0
                
        if self.T_inter[n,2]<0 and self.T_inter[n,4]>0:
            
            p1=0
            p2=1
            p3=1
            p4=0
            
                            
            if self.run_z1==True:
          
              #  print("cas 1")
                z1=-self.T_inter[n,2]*e_c/(self.Tl[n,0]-self.T_inter[n,2])
                z2=-self.T_inter[n,2]*e_c/(self.Tl[n,1]-self.T_inter[n,2])
                
                self.zfu[n,0]=Z1_2+(self.rho_l/self.rho_s)*(self.Cp_l/L)*self.T_inter[n,2]*z1/2
                self.zfu[n,1]=Z1_2-(self.rho_l/self.rho_s)*(self.Cp_l/L)*self.T_inter[n,2]*z2/2
      
        if self.T_inter[n,2]>0 and self.T_inter[n,4]<0  : 
                  
            p1=1
            p2=0
            p3=0
            p4=1
        
                
            if self.run_z3==True:
      
                z3=-self.T_inter[n,4]*e_c/(self.Tl[n,2]-self.T_inter[n,4])
                z4=-self.T_inter[n,4]*e_c/(self.Tl[n,3]-self.T_inter[n,4])
    
                
                self.zfu[n,2]=Z3_4+(self.rho_l/self.rho_s)*(self.Cp_l/L)*self.T_inter[n,4]*z3/2
                self.zfu[n,3]=Z3_4-(self.rho_l/self.rho_s)*(self.Cp_l/L)*self.T_inter[n,4]*z4/2
                    
            
        if self.T_inter[n,2]<0 and self.T_inter[n,4]<0  : 
        
            p1=0
            p2=1
            p3=0
            p4=1
          
            if self.run_z1==True:

                z1=-self.T_inter[n,2]*e_c/(self.Tl[n,0]-self.T_inter[n,2])
                z2=-self.T_inter[n,2]*e_c/(self.Tl[n,1]-self.T_inter[n,2])
                
                self.zfu[n,0]=Z1_2+(self.rho_l/self.rho_s)*(self.Cp_l/L)*self.T_inter[n,2]*z1/2
                self.zfu[n,1]=Z1_2-(self.rho_l/self.rho_s)*(self.Cp_l/L)*self.T_inter[n,2]*z2/2
            if self.run_z3==True:
                z3=-self.T_inter[n,4]*e_c/(self.Tl[n,2]-self.T_inter[n,4])
                z4=-self.T_inter[n,4]*e_c/(self.Tl[n,3]-self.T_inter[n,4])
    
                
                self.zfu[n,2]=Z3_4+(self.rho_l/self.rho_s)*(self.Cp_l/L)*self.T_inter[n,4]*z3/2
                self.zfu[n,3]=Z3_4-(self.rho_l/self.rho_s)*(self.Cp_l/L)*self.T_inter[n,4]*z4/2
                

          
        if dro==1:
            
            
            dromo.iterate(n+1,self.Tinj_dro[n]+kelvin,qdro_u)
            
            self.Tsor_dro[n+1]=dromo.T[n+1,1,-1]-kelvin
            
            if self.Tsor_dro[n+1]< self.T_inter[n,3] or self.Tsor_dro[n+1]< self.T_inter[n,5]:
                
                self.agenda_dro[n+1]=0
                initial=(self.T_be[n],self.T_inter[n,0],self.T_iso1[n],self.Tl[n,1],self.Ts[n,0]\
                 ,self.Tl[n,0],self.zfu[n,0],self.T_inter[n,2],self.Ts[n,1],self.Tl[n,1]\
                ,self.zfu[n,1],self.T_inter[n,3],self.Ts[n,2],self.Tl[n,2],self.zfu[n,2]\
                ,self.T_inter[n,4],self.Ts[n,3],self.Tl[n,3],self.zfu[n,3]\
                ,self.T_inter[n,5],self.Tl[n,4],self.T_inter[n,6],self.T_iso2[n]\
                ,self.Tinj_pac[n,0],self.Tinj_pac[n,1],self.Tsor_pac[n,0],self.Tsor_pac[n,1]\
                ,self.Tinj_sto[n,0],self.Tinj_sto[n,1],self.Tsor_sto[n,0],self.Tsor_sto[n,1]) # conditions initilaes à chaque itération
     
                solution=fsolve(self.equations,initial,args=(n,p1,p2,p3,p4)) # Calcul de la solution du système à l'instant n+1
                            
                self.Tinj_dro[n+1] = self.Tsor_dro[n+1]
            else:
                    
                initial=(self.T_be[n],self.T_inter[n,0],self.T_iso1[n],self.Tl[n,1],self.Ts[n,0]\
                 ,self.Tl[n,0],self.zfu[n,0],self.T_inter[n,2],self.Ts[n,1],self.Tl[n,1]\
                ,self.zfu[n,1],self.T_inter[n,3],self.Ts[n,2],self.Tl[n,2],self.zfu[n,2]\
                ,self.T_inter[n,4],self.Ts[n,3],self.Tl[n,3],self.zfu[n,3]\
                ,self.T_inter[n,5],self.Tl[n,4],self.T_inter[n,6],self.T_iso2[n]\
                ,self.Tinj_pac[n,0],self.Tinj_pac[n,1],self.Tsor_pac[n,0],self.Tsor_pac[n,1]\
                ,self.Tinj_sto[n,0],self.Tinj_sto[n,1],self.Tsor_sto[n,0],self.Tsor_sto[n,1]) # conditions initilaes à chaque itération
                self.o=np.zeros(self.size)   
                solution=fsolve(self.equations,initial,args=(n,p1,p2,p3,p4)) # Calcul de la solution du système à l'instant n+1
    
                self.Tmoy_sor_sto[n+1]=(solution[29]+solution[30])/2    
                self.Tinj_dro[n+1] = self.Tsor_dro[n+1] - self.eff * (self.Tsor_dro[n+1] - self.Tmoy_sor_sto[n+1])
                
        else:
            
            dromo.iterate(n+1,self.Tinj_dro[n]+kelvin,qdro_u)
            self.Tsor_dro[n+1]=dromo.T[n+1,1,-1]-kelvin        
            self.agenda_dro[n+1]=0  
            initial=(self.T_be[n],self.T_inter[n,0],self.T_iso1[n],self.Tl[n,1],self.Ts[n,0]\
                 ,self.Tl[n,0],self.zfu[n,0],self.T_inter[n,2],self.Ts[n,1],self.Tl[n,1]\
                ,self.zfu[n,1],self.T_inter[n,3],self.Ts[n,2],self.Tl[n,2],self.zfu[n,2]\
                ,self.T_inter[n,4],self.Ts[n,3],self.Tl[n,3],self.zfu[n,3]\
                ,self.T_inter[n,5],self.Tl[n,4],self.T_inter[n,6],self.T_iso2[n]\
                ,self.Tinj_pac[n,0],self.Tinj_pac[n,1],self.Tsor_pac[n,0],self.Tsor_pac[n,1]\
                ,self.Tinj_sto[n,0],self.Tinj_sto[n,1],self.Tsor_sto[n,0],self.Tsor_sto[n,1]) # conditions initilaes à chaque itération
  
            solution=fsolve(self.equations,initial,args=(n,p1,p2,p3,p4)) # Calcul de la solution du système à l'instant n+1       
            self.Tinj_dro[n+1] = self.Tsor_dro[n+1]
      

        self.T_be[n+1]=solution[0]
        self.T_iso1[n+1]=solution[2]
        self.T_inter[n+1,:]=solution[[1,3,7,11,15,19,21]]
        self.Ts[n+1,0:4]=solution[[4,8,12,16]]
        self.Tl[n+1,:]=solution[[5,9,13,17,20]]
        self.zfu[n+1,0:4]=solution[[6,10,14,18]] 
        self.T_iso2[n+1]=solution[22]
        self.Tinj_pac[n+1,0]= solution[23]
        self.Tinj_pac[n+1,1]= solution[24]
        self.Tsor_pac[n+1,0]= solution[25]
        self.Tsor_pac[n+1,1]= solution[26]
        self.Tinj_sto[n+1,0]=solution[27]
        self.Tinj_sto[n+1,1]=solution[28]                          
        self.Tsor_sto[n+1,0]=solution[29]
        self.Tsor_sto[n+1,1]= solution[30]   
        
        self.Tmoy_sor_sto[n+1]=(self.Tsor_sto[n+1,0]+self.Tsor_sto[n+1,1])/2
        self.Tmoy_inj_sto[n+1]=(self.Tinj_sto[n+1,0]+self.Tinj_sto[n+1,1])/2
        #self.Tmoy_inj_sto[n+1] = self.Tmoy_sor_sto[n+1] + self.coeff * self.eff * (self.Tsor_dro[n+1] - self.Tmoy_sor_sto[n+1])
        self.Tmoy_sor_pac[n+1]=(self.Tsor_pac[n+1,0]+self.Tsor_pac[n+1,1])/2
        self.Tmoy_inj_pac[n+1]=(self.Tinj_pac[n+1,0]+self.Tinj_pac[n+1,1])/2
     
        if self.T_inter[n+1,2]<=0 and  self.T_inter[n,2]<=0: 
            
            self.run_z1=False
            
        else:
            
            self.run_z1=True
            
            
        if self.T_inter[n+1,4]<=0 and  self.T_inter[n,4]<=0: 
           
           self.run_z3=False
           
        else:
            
            self.run_z3=True
            
        #print(self.run_z1)
            
        self.P_1[n+1]=self.mpac*self.cpac*(self.Tinj_pac[n+1,0]-self.Tsor_pac[n+1,0])    
       # self.P_2[n+1]=(self.T_inter[n+1,2]-0.5*(self.Tinj_pac[n+1,0]+self.Tsor_pac[n+1,0]))/self.Rthf                   
        self.P_2[n+1]=self.mpac*self.cpac*(self.Tinj_pac[n+1,1]-self.Tsor_pac[n+1,1])  
    
    def Tsous_sol(self,z,i):
        """
        La température du sous-sol est une fonction sinusoidale, obtenue à partir de l'équation de la propagation de la chaleur dans le sous-sol.
        
        La pulsation w a été ajustée dans setPertes pour que la période soit une année
        
        z : profondeur en mètres
        
        i : indice temporel dans la discrétisation
            
        """
        a = self.ath_l


         
        self.tf=18*24*3600
        # w pulsation exprimée en s-1
        self.w = 2*math.pi / (8760*3600)
        Tmoy=11
        Tamp=9.1        
        self.Tmoy=Tmoy 
        
        self.Tamp=Tamp
        
        self.za = math.sqrt( 2 * a / self.w )
        
  
        factor = math.cos(self.w * (self.step * i - self.tf) + z/self.za ) * math.exp(z /self.za)
            
        return self.Tmoy - self.Tamp * factor