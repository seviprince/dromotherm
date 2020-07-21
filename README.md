![codecheck](https://github.com/seviprince/dromotherm/workflows/codecheck/badge.svg)
![couplage1](https://github.com/seviprince/dromotherm/workflows/couplage1/badge.svg)
![couplage2](https://github.com/seviprince/dromotherm/workflows/couplage2/badge.svg)
![documentation](https://github.com/seviprince/dromotherm/workflows/documentation/badge.svg)

# Codes et Notebooks liés au projet Dromotherm

http://www.dromotherm.com

## master Energie

Depuis mars 2020

Dimensionnement/modélisation/choix de l’implantation d’un démonstrateur sur le site du Bourget de Lac à proximité de Chambéry

## bibliothèque dromosense

Le répertoire `dromosense` est appelé à herberger la bibliothèque d'outils liés au projet dromotherm

```
git clone https://github.com/seviprince/dromotherm
cd dromotherm
```
Pour installer la bibliothèque en mode développement : `pip install -e ./` ou `python3 setup.py develop`

La documentation est gérée au fil de l'eau par un workflow spécifique : https://seviprince.github.io/dromotherm/

La documentation se construit automatiquement à partir du code à l'aide de l'outil [pdoc](https://pdoc3.github.io/pdoc/)
```
pdoc --html -o docs dromosense
```

pour installer pdoc :

```
pip install pdoc3
```

## Privilégier le travail en environnement virtuel

### Anaconda

Avec la distribution Anaconda, le travail se fera d'office en environnement virtuel

Toujours lancer la console depuis Anaconda

### Créer/exploiter manuellement son environnement virtuel

Création d'un environnement appelé dromo
```
cd ~
python3.7 -m venv dromo
```
Celà crée dans les dossiers de l'utilisateur un dossier `dromo` avec tous les fichiers liés à l'environnement virtuel

Pour activer l'environnement
```
source dromo/bin/activate
```

Celà devrait afficher `(dromo)` dans l'invité de commande

Pour installer les bibliothèques, une fois l'environnement activé :
```
pip install numpy scipy matplotlib ipython jupyter pandas sympy nose
```
pour installer RISE, outil permettant de transformer un notebook en présentation :
```
pip install RISE
```
pour lancer jupyter :
```
jupyter notebook 
```
### MACOS

Pour installer pip
```
sudo easy_install pip
```


