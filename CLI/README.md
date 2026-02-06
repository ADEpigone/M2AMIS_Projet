Le CLI

*  Récupérer (avec un get) de molécules en ligne
* Gestion, affichage des différentes molécules téléchargées / ça implique de gérer un petit répo, ou une petite DB SQlite !
* Test d’iso de 2 molécules (juste l'appel aux fonctions faites par la personne qui gère ça)
* Test de similarité de 2 molécules (pareil)
* Affichage de la molécule en terminal ?

Pour la lib du CLI par exemple on peut utiliser : 

- https://click.palletsprojects.com/en/stable/

Ou tout faire à la main


Un exemple de --help pourrait être :

```cmd
$ nom_prog --help
Usage: nom_prog [OPTIONS] COMMAND [ARGS]...

Commands:
  list    Liste molécules
  get     Récupère depuis ChEBI (lien ou code jsp)
  delete  Supprime une molécule
  graph   Export graphe
  iso     Test isomorphisme
  show    Affiche molécule
  ... les distances etc après on verra
```

Détailler les commandes, la liste etc :
