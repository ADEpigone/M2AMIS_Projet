Permet de calculer la similarité de deux molécules ce qui retourne une valeur entre 0 et 1.

Il existe plusieurs options pour calculer une fingerprint (https://schwallergroup.github.io/practical-programming-in-chemistry/tutorials/lecture_05/03_rdkit_fingerprints.html):

- Morgan (https://towardsdatascience.com/a-practical-introduction-to-the-use-of-molecular-fingerprints-in-drug-discovery-7f15021be2b1/ ): fonctionne sur le voisinnage d'un atome (regarde de manière circulaire) en fonctions du radius
- rdkit : fonctionne en énumérant des sous arbres dans la mol 

Les atomes vont être attribués un int hashé à partir de leur num atomique, nb plus proche voisins (non H), nb de liens attachés, masse atomique, nb d'H, si il est dans un cycle. Puis tout ces int vont être modif en fonction de l'algo qui va donné la fingerpint finale.

Et plusieurs distances :
- Tanimoto (indice de Jaccard) : utilisé en pratique pour les molécules, vu que les fingerprints représente la structure de la molécules, Tanimoto compare les points communs des deux fingerprints de plus vu qu'il y a beaucoup de 0 dedans (bcp de propriétés diff possible pour une mol) Tanimoto est bien car il s'en fou des 0.
- Dice similarity : The Dice index is the ratio of the bits in common to the arithmetic mean of the number of on bits in the two items.
- Cosine : The Cosine index is the ration of the bits in common to the geometric mean of the number of on bits in the two items.

Tanimoto > 0.84 = très probable qu'ils aient des propriétés similaires.

https://www.daylight.com/dayhtml/doc/theory/theory.finger.html#RTFToC88

distances dispo dans RDKit : Tanimoto, Dice, Cosine, Sokal, Russel, Kulczynski, McConnaughey, and Tversky.

