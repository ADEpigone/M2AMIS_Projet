Pistes :

*  Base de cycle/Longueur de cycles, une métrique basée sur ça
* Distance d’édition de graphe (c’est safe mais NP complet, on peut utiliser ça comme baseline dans des benchmark je pense typiquement)
* Non supervisé (Clustering sur des embeddings de graphes, pq pas ! [https://dl.acm.org/doi/epdf/10.1145/3711896.3737128](https://dl.acm.org/doi/epdf/10.1145/3711896.3737128 "https://dl.acm.org/doi/epdf/10.1145/3711896.3737128") , il y a aussi une implémentation python dispo. SEUL PB c’est jsp si ça marche bien sur coloration, mais ça marche bien sur les poids, donc on peut p-ê encoder des infos dedans)
* Du ML, on peut essayer de faire de la similarité sémantique, regarder des GNN à ce lvl ou même entraîner notre propre encodeur vers espace latent. MAIS ce serait similarité sémantique donc il faut des infos sur les propriétés des molécules en plus de propriétés du graphe de la molécule. On peut aussi rester sur du syntaxique mais je doute qu’on puisse avoir meilleur modèle que SIGEM ou jsp.
* Laplacien des graphes : [https://elasalle.github.io/research/slides/2021\_10\_07\_datashape\_seminar.pdf](https://elasalle.github.io/research/slides/2021_10_07_datashape_seminar.pdf "https://elasalle.github.io/research/slides/2021_10_07_datashape_seminar.pdf") on peut utiliser des distances un peu exotiques basées sur des propriétés mathématiques ! ça a l'air super
