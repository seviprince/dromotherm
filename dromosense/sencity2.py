import numpy as np
from dromosense.constantes import kelvin
import math
from scipy.optimize import fsolve
from donnees import  Rthf, Cp_f,S_b, Z0,Zbe_i,Zi_1,Z1_2,Z2_3,Z3_4,Z4_5,Z5_i,Zs_sol,Qgeo, Tf, e_c, rho_l,Cp_l,rho_s,L

class UndergroundStorage:
    
    def __init__(self,fname,size,nc,step):
        
        self.step=step
        self.nc=nc
        self.Tinj_sto=5*np.ones((size,2))
        self.Tsor_sto=5*np.ones((size,2))
        self.Tmoy_inj_sto=5*np.ones(size)
        self.Tmoy_sor_sto=5*np.ones(size)
        
        self.Tsor_pac=5*np.ones((size,2))
        self.Tinj_pac=5*np.ones((size,2))
        self.Tmoy_inj_sto=5*np.ones(size)
        self.Tmoy_sor_sto=5*np.ones(size)       
        self.Tinj_dro=5*np.ones(size)
        self.Tsor_dro=5*np.ones(size)  
        
        self.Ts=0*np.ones((size,nc))
        self.Tl=5*np.ones((size,nc))
        
        self.T_inter=5*np.ones((size,nc))
        self.T_be=5*np.ones(size)
        self.T_iso1=5*np.ones(size)
        self.T_iso2=5*np.ones(size)
        
        _input=np.loadtxt(fname)
        
        self.Text=_input[:,1]
        

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
        
        self.Pgeo=np.zeros(size)


    def set(self,eff,k,coeff,cste,lambda_be,ath_be,lambda_iso,ath_iso,lambda_l,ath_l,lambda_s,ath_s,msto,cpsto,mpac,cpac):
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
        self.ath_be=ath_be,self.lambda_iso=lambda_iso
        self.ath_iso=ath_iso
        self.lambda_l=lambda_l
        self.ath_l=ath_l
        self.lambda_s=lambda_s
        self.ath_s=ath_s

        self.cpsto=cpsto
  
        self.cpac=cpac
       
        self.mpac=mpac
        self.cpac=cpac
       
        
        """
        calcul de B et C - cf notebook
        """
        
        self.B = (msto * cpsto - k/2) * coeff * eff
        
        self.C = 1 - k / (2 * mpac *cpac)
        
        print("le k du système géothermique vaut {} W/K".format(k))
        print("coeff vaut {}".format(coeff))
        print("B vaut {} W/K".format(self.B))
        print("C vaut {}".format(self.C))        

        
    
    def coeffs(self,T,Tsup,Tinf,Zsup,Zinf):
        
        """
        La température de chaque couche est sous une forme polynomiale T=az^2+b*z+c
        Cette fonction  calcule les coefficients a, b et c
        """
        
        """
    
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
        
        
        
        if self.Zinf==self.Zsup:
        
            return 0,0,0
    
        if self.Zinf==0.0 or self.Zsup==0.0:
            self.a=3*(-2*self.T+self.Tsup+self.Tinf)/(self.Zsup-self.Zinf)**2
            self.b=2*(3*self.T-self.self.Tsup-2*self.Tinf)/(self.Zsup-self.Zinf)
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
    
            self.c=(Y2*Z3-Y3*Z2)/(Y2*W3-Y3*W2)
            self.b=(Z2-W2*self.c)/Y2
            self.a=(Tinf-y1*self.b-self.c)/x1
            # On peut aussi utiliser la fonction linalg de numpy mais elle est coûteuse en temps d'éxécution     
        return self.a,self.b,self.c
    
    
    def T_interface(self,Tsup1,T1,Tinter,T2,Tinf2,Zsup1,Zinter,Zinf2,lambda1,lambda2,**kwargs):
   
        """
        
    
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
    
        Returns
        -------
        None.
    
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
                
                A=(1/(2*self.m*Cp_f)+Rthf)*S_b
                B=-self.lambda1*(2*self.a1*self.Zinter+self.b1)+self.lambda2*(2*self.a2*self.Zinter+self.b2)
                return self.Tinter-self.Tin-A*B


    def equations(self,Y,*p):
            
        
       """
        Définit le système des  équations décrivant les transferts thermiques dans le stockage avec ou
        
        sans changement de phase
    
        Y: vecteur contenant les 23 inconnues à l'instant n+1 (Tbe,Tbe_i,Tiso1,Ts1,Tl1,zfu1,....,Tiso2)
        
        Tout_pac=Tinj_pac
        Tin_pac=Tsor_pac
        
    
        """ 
        
       (Tbe,Tbe_i,Tiso1,Ti_1,Ts1,Tl1,zfu1,T1_2,Ts2,Tl2,zfu2,T2_3,Ts3,Tl3,zfu3,T3_4,Ts4,Tl4,zfu4,T4_5,Tl5,T5_i,Tiso2,Tinj_pac1,Tinj_pac2,Tsor_pac1,Tsor_pac2,Tsor_sto1,Tsor_sto2)=Y
       
       (n,p1,p2,p3,p4)=p
       a_be,b_be,c_be=self.coeffs(Tbe,self.Text,Tbe_i,Z0,Zbe_i)
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
       a_iso2,b_iso2,c_iso2=self.coeffs(Tiso2,T5_i,self.Tsous_sol(Zs_sol,n+1),Z5_i,Zs_sol) 
         
         
       #équations
         
       self.eq_be=Tbe-self.T_be[n]-2*self.ath_be*a_be*self.step
       self.eq_bi=self.T_interface(self.Text[n+1],Tbe,Tbe_i,Tiso1,Ti_1,Z0,Zbe_i,Zi_1,self.lambda_be,self.lambda_iso)
       self.eq_iso1=Tiso1-self.T_iso1[n]-2*self.ath_iso*a_iso1*self.step
         
       self.eqi_1=self.T_interface(Tbe_i,Tiso1,Ti_1,Tl1,T1_2*p1+p2*Tf,Zbe_i,Zi_1,Z1_2*p1+p2*zfu1,self.lambda_iso,self.lambda_l)
         
       self.eq_l1=Tl1-self.Tl[n,0]-2*self.ath_l*a_1*self.step
       self.eq_s1=(Ts1-self.Ts[n,0]-2*self.ath_s*a_p1*self.step)*p2
       self.eq_zfu1=(zfu1-self.zfu[n,0]-self.cste*(self.lambda_s*(2*a_p1*zfu1+b_p1)-self.lambda_l*(2*a_1*zfu1+b_1)))*p2
         
       self.eq_12=self.T_interface(Ti_1*p1+p2*Tf,Tl1*p1+p2*Ts1, T1_2,p1*Tl2+p2*Ts2, p1*T2_3+p2*Tf, p1*Zi_1+p2*zfu1, Z1_2, p1*Z2_3+p2*zfu2,p1*self.lambda_l+p2*self.lambda_s,p1*self.lambda_l+p2*self.lambda_s,m=self.mpac*self.agenda_pac[n+1],Tin=self.Tsor_pac[n+1,0])
            
       self.eq_l2=Tl2-self.Tl[n,1]-2*self.ath_l*a_2*self.step
       self.eq_s2=(Ts2-self.Ts[n,1]-2*self.ath_s*a_p2*self.step)*p2
       self.eq_zfu2=(zfu2-self.zfu[n,1]-self.cste*(self.lambda_s*(2*a_p2*zfu2+b_p2)-self.lambda_l*(2*a_2*zfu2+b_2)))*p2  
           
       self.eq_23=self.T_interface(p1*T1_2+p2*Tf, Tl2, T2_3, Tl3, p3*T3_4+p4*Tf, p1*Z1_2+p2*zfu2, Z2_3, p3*Z3_4+p4*zfu3, self.lambda_l, self.lambda_l,m=self.msto*self.agenda_dro[n+1],Tin=self.Tinj_sto[n+1,0])
         
       self.eq_l3=Tl3-self.Tl[n,2]-2*self.ath_l*a_3*self.step
       self.eq_s3=(Ts3-self.Ts[n,2]-2*self.ath_s*a_p3*self.step)*p4
       self.eq_zfu3=(zfu3-self.zfu[n,2]-self.cste*(self.lambda_s*(2*a_p3*zfu3+b_p3)-self.lambda_l*(2*a_3*zfu3+b_3)))*p4
         
       self.eq_34=self.T_interface(p3*T2_3+p4*Tf, p3*Tl3+p4*Ts3, T3_4,p3*Tl4+p4*Ts4, p3*T4_5+p4*Tf, p3*Z2_3+p4*zfu3, Z3_4, p3*Z4_5+p4*zfu4, p3*self.lambda_l+p4*self.lambda_s,  p3*self.lambda_l+p4*self.lambda_s,m=self.mpac*self.agenda_pac[n+1],Tin=self.Tsor_pac[n+1,1])
         
       self.eq_l4=Tl4-self.Tl[n,3]-2*self.ath_l*a_4*self.step
       self.eq_s4=(Ts4-self.Ts[n,3]-2*self.ath_s*a_p4*self.step)*p4
       self.eq_zfu4=(zfu4-self.zfu[n,3]-self.cste*(self.lambda_s*(2*a_p4*zfu4+b_p4)-self.lambda_l*(2*a_4*zfu4+b_4)))*p4   
         
       self.eq_45=self.T_interface(p3*T3_4+p4*Tf, Tl4, T4_5, Tl5, T5_i, p3*Z3_4+p4*zfu4, Z4_5, Z5_i, self.lambda_l, self.lambda_l,m=self.msto*self.agenda_dro[n+1],Tin=self.Tinj_sto[n+1,1])
         
       self.eq_l5=Tl5-self.Tl[n,4]-2*self.ath_l*a_5*self.step

         
       self.eq_5i=self.T_interface(T4_5,Tl5,T5_i,Tiso2,self.Tsous_sol(Zs_sol,n+1),Z4_5,Z5_i,Zs_sol,self.lambda_l,self.lambda_iso)
         
       # isolant 2
         
       self.eq_iso2=Tiso2-self.T_iso2[n]-2*self.ath_iso*a_iso2*self.step
       
       # Les fluides
       A=self.mpac*self.agenda_pac[n+1]*Rthf-0.5
       B=self.mpac*self.agenda_pac[n+1]*Rthf+0.5  
       A_prime=self.mpac*self.agenda_pac[n+1]*Rthf-0.5
       B_prime=self.mpac*self.agenda_pac[n+1]*Rthf+0.5            
       self.eq_inj_pac1=Tinj_pac1-(T1_2+ A*Tsor_pac1)/B  
       self.eq_inj_pac2=Tinj_pac2-(T3_4+ A*Tsor_pac2)/B   
       self.eq_sor_pac1=Tsor_pac1-Tinj_pac1+self.Pgeo/(2*self.mpac*self.agenda_pac[n+1]*Cp_f)
       self.eq_sor_pac2=Tsor_pac2-Tinj_pac2+self.Pgeo/(2*self.mpac*self.agenda_pac[n+1]*Cp_f)
          
       self.eq_sor_sto1=Tsor_sto1-(T2_3+ A_prime*self.Tin_dro1)/B_prime  
       self.eq_sor_sto2=Tsor_sto2-(T4_5+ A_prime*self.Tin_dro2)/B_prime       
       
    
       return self.eq_be,self.eq_bi,self.eq_iso1,self.eqi_1,self.eq_s1,self.eq_l1,self.eq_zfu1,self.eq_12,self.eq_s2,self.eq_l2,self.eq_zfu2,self.eq_23,self.eq_s3,self.eq_l3,self.eq_zfu3,self.eq_34,self.self.eq_s4,self.eq_l4,self.eq_zfu4,self.eq_45,self.eq_l5,self.eq_5i,self.eq_iso2, self.eq_inj_pac1, self.eq_inj_pac2, self.eq_sor_pac1, self.eq_sor_pac2,self.eq_sor_sto1,self.eq_sor_sto2
    #----------------------------------------------------------------------------------------------------------------------------------
      








    def StockLoop(self,i):
        """
        réalise une itération sur la température du stockage
        
        calcule la dérivée de la température du massif de stockage en K/s ou °C/s
        
        retourne la valeur de Tsable[i+1]
        
        4 cas distincts :
        
        1) appel d'énergie en provenance du bâtiment + dromotherme en marche
        
        2) appel d'énergie en provenance du bâtiment + dromotherme arrêté
        
        3) pas d'appel d'énergie en provenance du bâtiment + dromotherme en marche
        
        4) pas d'appel d'énergie en provenance du bâtiment + dromotherme arrêté
        
        modélise les pertes dans le stockage si la méthode setPertes a été utilisée
        """
        
        
        
        
        if self.simPertes==0:
            
            self.pertes[i]=0
        
        
        else:
    
            pertes_Laterales=self.SL_iso*(self.Tsable[i]-self.Tsous_sol(-1.8/2,i))
            
            pertes_Base=self.SB_iso*(self.Tsable[i]-self.Tsous_sol(-1.8,i))
            
            self.pertes[i]=self.u_th*(pertes_Laterales+ pertes_Base)
        
        pac=self.agenda_pac[i]
        dro=self.agenda_dro[i]
        if pac==1 and dro==1:
            der = (self.msto * self.cpsto * (self.Tinj_sto[i] - self.Tsor_sto[i]) - self.Pgeo[i]-self.pertes[i]) / (self.msable * self.cpsable)
        if pac==1 and dro==0:
            der = (- self.Pgeo[i]-self.pertes[i])/ (self.msable * self.cpsable)
        if pac==0 and dro==1:
            der = (self.msto * self.cpsto * (self.Tinj_sto[i] - self.Tsor_sto[i])-self.pertes[i]) / (self.msable * self.cpsable)
        if pac==0 and dro==0:
            der = -self.pertes[i]/(self.msable * self.cpsable)
            #der =0
        
        self.diff[i+1]=der
        
        ## schéma de discrétisation
        return self.Tsable[i]+self.step*der

    def SystemLoop(self,i,qdro_u,dromo):
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
        n=i
        
        if self.T_inter[n,2]>0 and self.T_inter[n,4]>0:
            
            p1=1
            p2=0
            p3=1
            p4=0
                
        if self.T_inter[n,2]<0 and self.T_inter[n,4]>0:
            
            p1=1
            p2=0
            p3=0
            p4=1
                                
            if self.run_z1==True:
          
       
                z1=-self.T_inter[n,2]*e_c/(self.Tl[n,0]-self.T_inter[n,2])
                z2=-self.T_inter[n,2]*e_c/(self.Tl[n,1]-self.T_inter[n,2])
                
                self.zfu[n,0]=Z1_2+(rho_l/rho_s)*(Cp_l/L)*self.T_inter[n,2]*z1/2
                self.zfu[n,1]=Z1_2-(rho_l/rho_s)*(Cp_l/L)*self.T_inter[n,2]*z2/2
           
        
           
            if self.T_inter[n,2]>0 and self.T_inter[n,4]<0  : 
                  
                p1=1
                p2=0
                p3=0
                p4=1
        
                
                if self.run_z3==True:
         
                    z3=-self.T_inter[n,4]*e_c/(self.Tl[n,2]-self.T_inter[n,4])
                    z4=-self.T_inter[n,4]*e_c/(self.Tl[n,3]-self.T_inter[n,4])
        
                    
                    self.zfu[n,2]=Z3_4+(rho_l/rho_s)*(Cp_l/L)*self.T_inter[n,4]*z3/2
                    self.zfu[n,3]=Z3_4-(rho_l/rho_s)*(Cp_l/L)*self.T_inter[n,4]*z4/2
                    
            
            if self.T_inter[n,2]<0 and self.T_inter[n,4]<0  : 
            
                p1=0
                p2=1
                p3=0
                p4=1
               
                if self.run_z1==True:

                    z1=-self.T_inter[n,2]*e_c/(self.Tl[n,0]-self.T_inter[n,2])
                    z2=-self.T_inter[n,2]*e_c/(self.Tl[n,1]-self.T_inter[n,2])
                    
                    self.zfu[n,0]=Z1_2+(rho_l/rho_s)*(Cp_l/L)*self.T_inter[n,2]*z1/2
                    self.zfu[n,1]=Z1_2-(rho_l/rho_s)*(Cp_l/L)*self.T_inter[n,2]*z2/2
                if self.run_z3==True:
                    z3=-self.T_inter[n,4]*e_c/(self.Tl[n,2]-self.T_inter[n,4])
                    z4=-self.T_inter[n,4]*e_c/(self.Tl[n,3]-self.T_inter[n,4])
        
                    
                    self.zfu[n,2]=Z3_4+(rho_l/rho_s)*(Cp_l/L)*self.T_inter[n,4]*z3/2
                    self.zfu[n,3]=Z3_4-(rho_l/rho_s)*(Cp_l/L)*self.T_inter[n,4]*z4/2
                       
        initial=(self.T_be[n],self.T_inter[n,0],self.T_iso1[n],self.Tl[n,1],self.Ts[n,0],self.Tl[n,0],self.zfu[n,0],self.T_inter[n,2],self.Ts[n,1],self.Tl[n,1],self.zfu[n,1],self.T_inter[n,3],self.Ts[n,2],self.Tl[n,2],self.zfu[n,2],self.T_inter[n,4],self.Ts[n,3],self.Tl[n,3],self.zfu[n,3],self.T_inter[n,5],self.Tl[n,4],self.T_inter[n,6],self.T_iso2[n],self.Tinj_pac[n,0],self.Tinj_pac[n,1],self.Tsor_pac[n,0],self.Tsor_pac[n,1],self.Tsor_sto[n,0],self.Tinj_sto[n,1]) # conditions initilaes à chaque itération

        solution=fsolve(self.equations,initial,args=(n,p1,p2,p3,p4)) # Calcul de la solution du système à l'instant n+1
   
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
        self.Tsor_sto[n+1,0]=solution[27]
        self.Tsor_sto[n+1,1]=solution[28]                          
            
        
        # étape 2        
        dro=self.agenda_dro[n]
        pac=self.agenda_pac[n]
        
        y = np.mean(self.Tl[n+1,:])
    

        if pac==0:
            self.Tinj_pac[n+1,0]= self.T_inter[n+1,2]
            self.Tinj_pac[n+1,1]= self.T_inter[n+1,4]
            self.Tsor_pac[n+1,0]= self.T_inter[n+1,2]
            self.Tsor_pac[n+1,1]= self.T_inter[n+1,4]    
            self.Tmoy_sor_sto[n+1]=(self.Tsor_sto[n+1,0]+self.Tsor_sto[n+1,1])/2 
            self.Tmoy_inj_sto[n+1]=(self.Tinj_sto[n+1,0]+self.Tinj_sto[n+1,1])/2
        else:
            self.Tinj_pac[n+1,0]= solution[23]
            self.Tinj_pac[n+1,1]= solution[24]
            self.Tsor_pac[n+1,0]= solution[25]
            self.Tsor_pac[n+1,1]= solution[26]            
            self.Tmoy_sor_sto[n+1]=(self.Tsor_sto[n+1,0]+self.Tsor_sto[n+1,1])/2 
            self.Tmoy_inj_sto[n+1]=(self.Tinj_sto[n+1,0]+self.Tinj_sto[n+1,1])/2               
        
        # étape 3
        if dro == 1:
            dromo.iterate(n+1,self.Tinj_dro[n]+kelvin,qdro_u)
            self.Tsor_dro[n+1]=dromo.T[n+1,1,-1]-kelvin
            if self.Tsor_dro[n+1] < y :
                #print("step {} y vaut {} et prev vaut {}".format(i,y,Tsor_dro[i]))
                self.agenda_dro[n+1] = 0
                self.Tinj_dro[n+1] = self.Tsor_dro[n+1]
                self.Tinj_sto[n+1,0] = self.T_inter[n+1,3]
                self.Tinj_sto[n+1,1] = self.T_inter[n+1,5]
                self.Tsor_sto[n+1,0] = self.T_inter[n+1,3]
                self.Tsor_sto[n+1,1] = self.T_inter[n+1,5]
            else :
                
                self.Tmoy_sor_sto[n+1]=(self.Tsor_sto[n+1,0]+self.Tsor_sto[n+1,1])/2
                self.Tmoy_inj_sto[n+1] = self.Tmoy_sor_sto[n+1] + self.coeff * self.eff * (self.Tsor_dro[n+1] - self.Tmoy_sor_sto[n+1])
                self.Tinj_dro[n+1] = self.Tmoy_sor_dro[n+1] - self.eff * (self.Tsor_dro[n+1] - self.Tmoy_sor_sto[n+1])
            
        else:
            dromo.iterate(n+1,self.Tinj_dro[n]+kelvin,0)
            self.Tinj_dro[n+1] = self.Tsor_dro[n+1]
            self.Tinj_sto[n+1,0] = self.T_inter[n+1,3]
            self.Tinj_sto[n+1,1] = self.T_inter[n+1,5]
            self.Tsor_sto[n+1,0] = self.T_inter[n+1,3]
            self.Tsor_sto[n+1,1] = self.T_inter[n+1,5]

        
             
    def setPertes(self,Tmoy,Tamp,lambda_sable,rho_sable,e_iso,SL_iso,SB_iso,lambda_iso):   
        """
        paramètres à prendre en compte pour la modélisation des pertes
         
        Tmoy : Température moyenne annuelle du sous-sol en °C
        
        Tamp : l'amplitude annuelle de la température i.e Tmax-Tmin en °C
        
        tf : jour correspondant à la température minimale du sous-sol (18 janvier à Chambéry) en secondes écoulées depuis le début de l'année
        
        lambda_sable : conductivité thermique du sous-sol en W/(K.m)
        rho_sable : masse volumique du sable en kg/m^3
        
        lambda_iso : conductivité thermique de l'isolant autour du stockage en W/(K.m)
        
        SL_iso : surface de l'isolant le long des parois latéraux du stockage en m2
        
        SB_iso : surface de l'isolant posée à la base du stockage en m2
        
        e_iso : épaisseur de l'isolant en m
        
        """
        
        self.simPertes=1
         
        self.tf=18*24*3600
        # w pulsation exprimée en s-1
        self.w = 2*math.pi / (8760*3600)
         
        self.Tmoy=Tmoy 
        
        self.Tamp=Tamp
        
        self.SL_iso=SL_iso
        
        self.SB_iso=SB_iso

        # a diffusivité thermique du sous-sol, en m2/s
        a = lambda_sable / (self.cpsable * rho_sable)

        self.za = math.sqrt( 2 * a / self.w )
        
        self.u_th = lambda_iso / e_iso
        
    
    def Tsous_sol(self,z,i):
        """
        La température du sous-sol est une fonction sinusoidale, obtenue à partir de l'équation de la propagation de la chaleur dans le sous-sol.
        
        La pulsation w a été ajustée dans setPertes pour que la période soit une année
        
        z : profondeur en mètres
        
        i : indice temporel dans la discrétisation
            
        """
        
        factor = math.cos(self.w * (self.step * i - self.tf) + z/self.za ) * math.exp(z /self.za)
            
        return self.Tmoy - self.Tamp * factor