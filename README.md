## Sujet: Cartes Topologiques de Kohonen et coudes protéiques

Nicolas Silva (silva.nicolas.j@gmail.com)<br/>
Université Paris Diderot - 2019/2020 - Projet Long

__Objectif__

Le but de ce projet est d'implémenter un modèle capable de classifier des données basées sur des angles et amener un apprentissage différent dans
les cartes de Kohonen (SOM ou Self-Organizing Maps).<br/><br/>

Exécution sous l'environnement python3

#### Répertoires et localisation des fichiers
*********************************************

Le répertoire contient :<br/>
	- un dossier *data* contenant le fichier de données *.db* utilisé pour les analyses.<br/>
	- un dossier *doc* contenant le rapport du projet.<br/>
	- un dossier *results* contenant les différents résultats obtenus (fréquences des acides aminées, des structures, cartes de Kohonen, ...).<br/>
	- un dossier *src* contenant les scripts utilisés.<br/>
	- un fichier *.yml* permettant de recréer l'environnement conda utilisé.<br/>

#### Programme
**************

Il y a un progamme à exécuter : le programme *kohonen_BetaTurn_python3.py* qui prend en argument le nom du fichier de données à analyser.<br/>
	
'''Ex: *python3 kohonen_BetaTurn_python3.py -f shcullpdb_pc20_res2.5_R1.0_d050827_chains2722.db*'''

#### Résultats
*************************************

Ils peuvent être trouvés dans le dossier *results*.<br/>
Sont présentés les différentes fréquences d'acides aminés, de structures secondaires et de type de coudes Beta, ainsi que les résultats de l'apprentissage via les cartes de Kohonen.<br/>
La partie des résultats est décrite dans le rapport dans le dossier *doc*.

#### Modules python importés
*****************************

Os, Argparse, Math, Matplotlib.pyplot, Numpy, Random
