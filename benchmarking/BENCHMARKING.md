# Benchmarking – Architecture et état actuel

## Objectif

Ce module implémente une **première architecture de benchmarking** pour comparer différentes distances entre molécules, en cohérence avec l’approche théorique présentée (kNN, monotonie, séparation).

À ce stade, il s’agit volontairement d’un **squelette fonctionnel**, dont l’objectif principal est de :
- fixer une **architecture propre et extensible**,
- découpler les **distances** des **méthodes de benchmarking**,
- permettre l’ajout futur de distances et de datasets réels **sans modifier la structure du code**.

---

## Interface `Distance`

Une interface abstraite `Distance` est définie avec une méthode centrale :

```python
calculer_distance(mol1, mol2) -> float
```

- Cette méthode retourne une **valeur réelle** représentant une dissimilarité
  (plus la valeur est faible, plus les molécules sont similaires).
- Le type des molécules (`mol1`, `mol2`) est volontairement générique à ce stade.
- Dans la version finale du projet, ces objets correspondront à des **molécules parsées sous forme de graphes**.

**Point clé :**  
Les méthodes de benchmarking ne dépendent que de cette interface, et non de l’implémentation concrète de la distance.

---

## Distances actuelles

Une implémentation factice de la distance est fournie à des fins de test :
- elle retourne une valeur pseudo-aléatoire,
- elle permet de valider le pipeline de benchmarking indépendamment des distances finales.

Cette implémentation sera remplacée ultérieurement par des distances concrètes
sans impact sur les méthodes de benchmarking.

---

## Factory Method

La création des distances est centralisée via un mécanisme de **factory method** :
- la distance est sélectionnée à partir de son nom,
- le code de benchmarking reste indépendant des implémentations concrètes.

---

## Méthodes de benchmarking

Deux méthodes de benchmarking sont actuellement implémentées.

### 1. Benchmark k-Nearest Neighbors (kNN)

Cette méthode :
- calcule les k plus proches voisins d’une molécule,
- compare le voisinage obtenu avec un voisinage attendu lorsque la vérité terrain est disponible,
- permet le calcul de métriques telles que :
  - précision@k,
  - rappel@k (lorsqu’il est défini),
  - et, à terme, des métriques de stabilité.

Les datasets utilisés sont pour l’instant **simples et artificiels**,
et servent uniquement à tester la logique de benchmarking.

---

### 2. Benchmark Monotonie / Séparation

Cette méthode teste des **propriétés qualitatives attendues** sur des jeux de données de type :

```
[seed, e1, e2, e3, random]
```

- **Monotonie** : la distance à la seed doit augmenter avec le nombre de modifications.
- **Séparation** : les variantes doivent être plus proches de la seed que des molécules non pertinentes.

Cette approche est conceptuellement proche d’une **triplet loss**, mais utilisée ici uniquement comme outil d’évaluation.

---

## Datasets : état actuel et évolution prévue

À l’heure actuelle :
- les datasets sont **artificiels et figés**,
- ils servent uniquement à valider le fonctionnement du pipeline.

À terme :
- des datasets concrets et réalistes seront intégrés,
- adaptés spécifiquement à chaque méthode de benchmarking,
- sans modification de l’architecture ni des interfaces existantes.

cette partie dépend fortement des dataset de benchmarking, une fois que j'ai les deux dataset nécéssaires pour les deux méthodes les choses vont être 100% concrètes 
