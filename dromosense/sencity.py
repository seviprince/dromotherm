import numpy as np
from dromosense.constantes import kelvin


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
    
    agenda_dro et agenda_pac : agendas de fonctionnement du dromotherme et de la PAC
    
    Pgeo : puissance géothermique à développer pour satisfaire le besoin du bâtiment en W
    
    par exemple, avec une PAC de COP 3, la puissance géothermique à développer vaut 2/3 du besoin total du bâtiment
    """
    
    def __init__(self,size,step):
        """
        size : nombre de points pour la discrétisation temporelle
        
        step : pas de temps en secondes
        """
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
        """
        pac=self.agenda_pac[i]
        dro=self.agenda_dro[i]
        if pac==1 and dro==1:
            der = (self.msto * self.cpsto * (self.Tinj_sto[i] - self.Tsor_sto[i]) - self.Pgeo[i]) / (self.msable * self.cpsable)
        if pac==1 and dro==0:
            der = - self.Pgeo[i] / (self.msable * self.cpsable)
        if pac==0 and dro==1:
            der = self.msto * self.cpsto * (self.Tinj_sto[i] - self.Tsor_sto[i]) / (self.msable * self.cpsable)
        if pac==0 and dro==0:
            der = 0
        
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
        
        ```
        cas 1: le dromotherme est en marche 
        le fluide circule avec un débit unitaire qdro_u
        Test = Tsor_dro > Tsable ?
          
        Test négatif : pas d'échange d'énergie entre la route et le stock, 
        cf dromotherme à l'arrêt + on passe la valeur de agenda_dro[i] à 0
          
        Test positif : alimentation du stockage
        ```

        ```        
        cas 2: le dromotherme est à l'arrêt
        le débit est nul
        l'échangeur de séparation de réseau ne tourne pas
        
        pas de prélèvement par l'échangeur de séparation de réseau
        Tinj_dro[i] = Tsor_dro[i]
           
        fonctionnement à perte nulle pour le stockage
        Tsor_sto[i]=Tsor_sto[i-1] et Tinj_sto[i]=Tinj_sto[i-1]
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

