Permet de calculer la similarité de deux molécules ce qui retourne une valeur entre 0 et 1.

Il existe plusieurs options pour calculer une fingerprint :
- Morgan :
- rdkit :

Et plusieurs distances :
- Tanimoto (indice de Jaccard) : utilisé en pratique pour les molécules, je crois.
- Dice similarity : je sais pas encore pk en particulier.
- Cosine : pareil.

distances dispo dans RDKit : Tanimoto, Dice, Cosine, Sokal, Russel, Kulczynski, McConnaughey, and Tversky.

The usual metric for similarity between atom-pair fingerprints is Dice similarity