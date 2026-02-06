# Un peu d'algèbre

Soit G un graphe, on note A sa matrice d'adjacence, on note L son laplacien.
Le laplacien est la matrice d'adjacence où l'on enlève le degré du sommet sur sa diagonale.

On peut alors regarder le spectre des graphes et en faire la distance, un exemple de distance serait alors :
    d(G,H) = || spectre(G) - spectre(H) || avec ||.|| une norme au choix (L1,L2 etc ... à voir ?)

Mais que faire avec les couleurs ? Pour les arêtes on peut juste mettre dans la matrice d'adjacence le nombre, comme c'est le nombre de liens.

Rapide : (pas tout capté)
https://people.mpi-inf.mpg.de/~mehlhorn/ftp/genWLpaper.pdf
Clair : 
https://arxiv.org/pdf/2201.07083


Graph kernel k(G,H) = distance entre les deux graphes.
On a k(G,H) = <phi(G), phi(H)> avec phi un embedder (souvent infini apparemment)
Il faut que l'espace embeddé soit un espace de Hilbert (a un produit scalaire et est complet sous la norme induite). On n'a aucune condition sur phi !

1-WL : 
Itérativement on va construire de nouveaux labels sur les graphes, pour chaque sommet on va faire hash(label sommet + liste labels voisins)
Ensuite pour calculer le kernel on regarde les différences
Avantages : 
- Complexité linéaire/quadratique (papier genWL)
- Simple
Désavantages :
- Je sais pas encore, moins bien que Deep Learning ?

Autres kernels:
Random Walk, les plus courts chemins entre labels
Avantages :
- Comparables à 1-WL ?
Désavantages :
- Lents

GNN :https://arxiv.org/pdf/1810.00826 -> GIN ~~ 1-WL 

k-WL : 1-WL mais on prend des tuples de taille k (et non 1)
n-WL : test d'isomorphisme

Heat Diffusion :
Ht = e^-tL (t +- un zoom, comme sigma dans RFF)

On peut essayer de vectoriser les labels et avoir Ht x Labels 
Et le kernel pourrait être la somme de ces vecteurs 
ON PEUT AUSSI : Essayer de faire apprendre la fonction, plutôt que de sommer

Avantage : 
- Jsp ça a l'air pas mal

Désavantages:
- Potentiellement super lent
