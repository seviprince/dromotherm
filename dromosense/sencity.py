import numpy as np
from dromosense.constantes import kelvin
import math

class SenCityOne:
    """
    Cette classe réalise un premier couplage du système complet et nourrit un ensemble de vecteurs numpy:
    
    Tinj_sto : température du fluide de charge à l'injection dans le stockage
    
    Tsor_sto : température du fluide de charge après transit dans le stockage
    
    Tsor_pac : température du fluide de décharge au sortir de la PAC à la réinjection dans le stockage
    
    Tinj_pac : température du fluide de décharge à l'injection dans la PAC
    
    Tsable : température du stockage/sable
    
    diff : valeur de la dérivée de la température du stockage en °C/s ou K/s
    
    Tinj_dro : température d'injection dans le dromotherme
    
    Tsor_dro : température de sortie du dromotherme
    
    On initialise Tinj_dro et Tsor_dro à 10
    
    pertes : pertes instantannées dans le stockage en W (prérequis : nécessite l'utilisation de setPertes)
    
    agenda_dro et agenda_pac : agendas de fonctionnement du dromotherme et de la PAC
    
    Pgeo : puissance géothermique à extraire pour satisfaire le besoin du bâtiment en W
    
    par exemple, avec une PAC de COP 3, la puissance géothermique à extraire vaut 2/3 du besoin total du bâtiment
    
    Pour l'utiliser :
    
    ```
    # instanciation
    # RSB = route stock bâtiment
    RSB=SenCityOne(meteo.shape[0],3600)
    # définition du système
    RSB.set(eff,k,coeff,msto,cpsto,msable,cpsable,mpac,cpac)
    # injection du besoin
    RSB.Pgeo = (COP-1) * besoin_total / COP
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
        RSB.SystemLoop(i,qdro_u,dromo)
    ```
    """
    
    def __init__(self,size,step):
        """
        size : nombre de points pour la discrétisation temporelle
        
        step : pas de temps en secondes
        """
        self.simPertes = 0
        self.step=step
        
        self.Tinj_sto=np.zeros(size)
        self.Tsor_sto=np.zeros(size)
        
        self.Tsor_pac=np.zeros(size)
        self.Tinj_pac=np.zeros(size)
        
        self.Tsable=np.zeros(size)
        self.diff=np.zeros(size)
        
        self.Tinj_dro=10*np.ones(size)
        self.Tsor_dro=10*np.ones(size)

        self.agenda_dro=np.zeros(size)
        self.agenda_pac=np.zeros(size)
        
        self.pertes=np.zeros(size)
        
        self.Pgeo=np.zeros(size)
        
    def set(self,eff,k,coeff,msto,cpsto,msable,cpsable,mpac,cpac):
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
        self.msto=msto
        self.cpsto=cpsto
        self.msable=msable
        self.cpsable=cpsable
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
    
            pertes_Laterales=self.SL_iso*(self.Tsable[i]-self.Tsous_sol(-1.125,i))
            
            pertes_Base=self.SB_iso*(self.Tsable[i]-self.Tsous_sol(-2.25,i))
            
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
        
        Test = Tsor_dro > Tsable ?
          
        Test négatif : goto cas 2 (dromotherme à l'arrêt) + on passe la valeur de agenda_dro[i] à 0
          
        Test positif : alimentation du stockage
        ```

        ********************************************************************************       
        cas 2: le dromotherme est à l'arrêt - débit nul
        
        l'échangeur de séparation de réseau ne tourne pas
        
        pas de prélèvement par l'échangeur de séparation de réseau
        
        `Tinj_dro[i] = Tsor_dro[i]`
        
        pas d'évolution des températures du stockage
        
        `Tsor_sto[i]=Tsor_sto[i-1]` et `Tinj_sto[i]=Tinj_sto[i-1]`
        ```
          
        """
        # étape 1    
        self.Tsable[i]=self.StockLoop(i-1)
        
        dro=self.agenda_dro[i]
        pac=self.agenda_pac[i]
        
        y = self.Tsable[i]
    
        # étape 2
        if pac == 1 :
            self.Tinj_pac[i] = y-self.C*self.Pgeo[i]/self.k
            self.Tsor_pac[i] = self.Tinj_pac[i]-self.Pgeo[i]/(self.mpac*self.cpac)
        else:
            self.Tinj_pac[i] = self.Tinj_pac[i-1]
            self.Tsor_pac[i] = self.Tsor_pac[i-1]
    
        # étape 3
        if dro == 1:
            dromo.iterate(i,self.Tinj_dro[i-1]+kelvin,qdro_u)
            self.Tsor_dro[i]=dromo.T[i,1,-1]-kelvin
            if self.Tsor_dro[i] < y :
                #print("step {} y vaut {} et prev vaut {}".format(i,y,Tsor_dro[i]))
                self.agenda_dro[i] = 0
                self.Tinj_dro[i] = self.Tsor_dro[i]
                self.Tinj_sto[i] = self.Tinj_sto[i-1] 
                self.Tsor_sto[i] = self.Tsor_sto[i-1]
            else :
                self.Tsor_sto[i] = ( self.k * y + self.B * self.Tsor_dro[i] ) / ( self.k + self.B)
                self.Tinj_sto[i] = self.Tsor_sto[i] + self.coeff * self.eff * (self.Tsor_dro[i] - self.Tsor_sto[i])
                self.Tinj_dro[i] = self.Tsor_dro[i] - self.eff * (self.Tsor_dro[i] - self.Tsor_sto[i])
            
        else:
            dromo.iterate(i,self.Tinj_dro[i-1]+kelvin,0)
            self.Tsor_dro[i] = dromo.T[i,1,-1]-kelvin
            self.Tinj_dro[i] = self.Tsor_dro[i]
            self.Tinj_sto[i] = self.Tinj_sto[i-1] 
            self.Tsor_sto[i] = self.Tsor_sto[i-1]

        
             
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
        La température du sous-sol est une fonction sinusoidale :
        - de période = une année
        - de pulsation w ,obtenue à partir de l'équation de la propagation de la chaleur dans le sous-sol.
        
        z : profondeur en mètres
        
        i : indice temporel dans la discrétisation
            
        """
        
        factor = math.cos(self.w * (self.step * i - self.tf) + z/self.za ) * math.exp(z /self.za)
            
        return self.Tmoy - self.Tamp * factor
    
    
    
    
    
        
    
    
    