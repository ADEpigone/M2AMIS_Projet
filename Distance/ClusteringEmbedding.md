
Ici l'idée est de prendre le graphe + données et de calculer un embedding dans un espace latent.
L'espace latent est un espace abstrait de "haute" dimension (128 à 512 ?).
Il faut alors définir une fonction f : G,D -> L.

Pour cela dans le cas non coloré on a : https://dl.acm.org/doi/epdf/10.1145/3711896.3737128

Sinon : https://github.com/benedekrozemberczki/graph2vec


Sinon fait main :

On peut entraîner un GNN pour embed : https://en.wikipedia.org/wiki/Graph_neural_network

Puis on peut entraîner un embedder qui prend embedd graphe + embed data et s'occupe d'embed les 2.
On entraîne l'embedder global, ce qui entraîne les 2 en même temps et garde la cohérence.

Avec la taille de nos molécules on peut sûrement avoir un modèle spécialisé avec ~ 300k params.

Il faut récupérer les morceaux importants pour les données, ils sont :
- jsp g pas fait


Avantages : 
- Syntaxique ET sémantique
- Potentiellement très rapide

Désavantages :
- Demande bcp de données, compliqué à train ?
- Stochastique
