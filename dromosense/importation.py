import numpy as np
import matplotlib.pyplot as plt
import os

verbose=False


def getCsvDatas(name, dir="../datas", preview=False):
    """
    IMPORTATION de données depuis un fichier CSV

    The timeserie must be a fixed interval timeserie

    première colonne : le temps en seconde

    colonnes 1 à x : données

    name: file name

    dir: relative folder's path of the file

    preview: boolean - do yu want to preview the datas as plot charts ?

    return :

    step : interval in seconds

    float_data : the matrix containing the datas
    """
    fpath = "{}/{}".format(dir,name)
    f = open (fpath)
    data = f.read()
    f.close()
    """
    1) remove blank lines which could pollute the end of the file
    2) split into a list of lines
    3) extract the header and remove it from the object lines
    """
    data = data.rstrip('\n\r')
    lines = data.split('\n')
    header = lines[0].split(',')
    lines = lines[1:]
    # evaluates the timestep in s
    step = float(lines[1].split(',')[0])-float(lines[0].split(',')[0])
    step=int(step)
    if verbose:
        print(header)
        print(len(lines))
        input("press any key")

    """
    float_data shape is (time,features)
    on crée un tenseur rempli de zéros:
    - la dimension 0 (ou premier axe du tenseur) va être le nombre de données
    - la dimension 1 (ou second axe du tenseur) va être len(header)-1
    En effet, on met de côté la colonne liée au temps, qui ne sert a rien vu qu'on travaille à intervalle de temps fixe
    Généralement, dans un tenseur, on dit que le premier axe (axis 0) est le "samples axis"
    Ensuite on parcourt l'objet lines, on splite au niveau des , et on ne conserve que la partie [1:], c'est à dire tous les éléments sauf la premier valeur qui est le temps
    la numérotation commence à 0, donc :
    - line.split(',')[0] est le temps
    - line.split(',')[1] est la température de l'air
    - line.split(',')[2] est la température du point de rosée
    etc, etc
    """
    float_data=np.zeros((len(lines),len(header)-1))
    for i,line in enumerate(lines):
        values = [float(x) for x in line.split(',')[1:]]
        float_data[i,:]=values

    print("data shape is {}".format(float_data.shape))
    print("interval is {}".format(step))

    """
    float_data est donc l'objet contenant toutes tes données sauf le temps....
    c'est la méthode idéale dans ton cas car toutes tes données sont dans un seul fichier
    effectivement on aurait pu faire différemment s'il y avait eu un seul fichier par type de données, ie un fichier pour la température de l'air, un pour la radiation....
    """
    if preview:
        xrange=np.arange(float_data.shape[0])
        ax1=plt.subplot(111)
        plt.title("test d'importation")
        plt.xlabel('Temps (en heure)')
        ax1.set_ylabel("°C")
        plt.plot(float_data[:,0], label="Text", color="blue")
        plt.legend(loc='upper left')
        ax2=ax1.twinx()
        ax2.set_ylabel("W/m2")
        plt.fill_between(xrange,float_data[:,5], label="ray. atm. (W/m2)", color="orange", alpha=0.4)
        plt.fill_between(xrange,float_data[:,4], label="glob. rad. (W/m2)", color="yellow", alpha=0.4)
        plt.legend(loc='upper right')
        plt.show()

    return step, float_data

def getTxtModelicaDatas(name, dir="../datas", preview=False):
    """
    IMPORTATION de données txt comme attendu par Modelica

    Attention, modelica permet souvent de travailler avec des fichiers à pas de temps variable

    la colonne temps est donc à conserver
    """
    fpath = "{}/{}".format(dir,name)
    f = open (fpath)
    data = f.read()
    f.close()
    data = data.rstrip('\n\r')
    lines = data.split('\n')
    header = lines[1].split('#Columns:')[1].rstrip('\t').split(';')
    lines = lines[2:]
    if verbose:
        print(header)
    line1=[float(x) for x in lines[0].split('\t')]

    float_data=np.zeros((len(lines),len(line1)))
    for i,line in enumerate(lines):
        values = [float(x) for x in line.split('\t')]
        float_data[i,:]=values

    print("data shape is {}".format(float_data.shape))

    if preview:
        ax1=plt.subplot(111)
        plt.title("test d'importation")
        plt.xlabel('secondes')
        ax1.set_ylabel("W/m2")
        for i in range(1,3):
            plt.plot(float_data[:,0],float_data[:,i], label=header[i].strip())
        plt.legend(loc='upper left')
        plt.show()

    return float_data
