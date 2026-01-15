# La distance d'édition d'un graphe

Cette distance est purement syntaxique, aucune sémantique.

On va regarder le nb d'éditions à faire pour arriver de l'un à l'autre, c'est leur distance dans l'espace des graphes.
On autorise :
- La suppression d'une arête
- La suppression d'un sommet (mais il doit être déconnecté)
- L'ajout d'une arête
- L'ajout d'un sommet


Pseudo-code :
```
def GED(G, H):

    best_cost = +inf

    #un etat = (G_cur, cout_deja_paye)
    initial_state = (G, 0)

    queue = { initial_state }

    while queue non vide:
        cur = queue.pop(0)

        #si on a déjà dépassé
        if cur.cout_deja_paye >= best_cost:
            continue

        if cur.G_cur est identique a H:
            best_cost = min(best_cost, cur.cout_deja_paye)
            continue
        #on récup quelque chose qui diffère
        diff = choisir_une_difference(cur.G_cur, H)

        successeurs = []

        if la diff est une arete en trop dans G_cur:
            successeurs.add(appliquer(supprimer arete diff, cur))

        if la diff est une arete manquante dans G_cur:
            successeurs.add(appliquer(ajouter arete diff, cur))

        if la diff est une sommet en trop dans G_cur:
            successeurs.add(appliquer(supprimer une arete incidente diff, cur))
            if sommet_degre_zero(diff.sommet, cur.G_cur):
                successeurs.add(appliquer(supprimer sommet diff, cur))

        if la diff est sommet manquant dans G_cur:
            successeurs.add(appliquer(ajouter sommet diff, cur))

        if la diff est une couleur de sommet :
            successeur.add(état avec la couleur changée) 
        if la diff est une couleur d'arête :
            successeur.add(état avec la couleur changée) 

        for s in successeur:
            if s.cout_deja_paye + lower_bound(s.G_cur, H) < best_cost:
                queue.add(s)

    return best_cost
```

Avantages :
- Déterministe et exact (syntaxiquement)
- Intuitif
- Bon baseline pour du benchmarking syntaxique

Désavantages :
- Purement syntaxique, aucune composante syntaxique
- Pb NP dur