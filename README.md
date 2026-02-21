# M2AMIS Projet

Dépôt du code pour le projet de l'UE Projet AMIS.

Ici réside le code utilisé dans le rapport du groupe B.

## Sommaire

1. Fonctionnalités Principales
2. Architecture du Projet
3. Installation
4. Démarrage Rapide (Quickstart)
5. Guide du CLI

---

## Fonctionnalités Principales

- **Récupération et cache local** : Téléchargement automatique depuis ChEBI et stockage dans une BDD SQLite.
- **Test d'isomorphisme** : Vérification d'isomorphisme exact entre deux molécules.
- **Comparaison de similarité** : Calcul de similarité entre graphes moléculaires via Cellular WL ou des fingerprints classiques (Morgan, RDKit).
- **Clustering et analyse** : Clustering hiérarchique, analyse de pureté des clusters et identification des familles dominantes.
- **Benchmark ESOL (`logS`)** : Évaluation des performances de prédiction de solubilité dans l'eau et de continuité du voisinage.
- **Corrélation Syntaxique/Sémantique** : Comparaison entre la similarité structurelle (fingerprints) et une similarité sémantique (ontologie ChEBI).

---

## Architecture du Projet

Le projet est structuré de manière modulaire pour séparer la logique métier, l'accès aux données et l'interface utilisateur (CLI).

```text
M2AMIS_Projet/
├── mol_cli.py                 # CLI à utiliser pour lancer le code
├── utils.py                   # Fonctions utilitaires
├── graph.py                   # Structure de données utilisée
│
├── Chebi/                     # Gestion des données
│   ├── CheBi2.py              # Interface SQLite avec backup API CheBi
│   └── ontology/              # Gestion des ontologies CheBi
│
├── cli_plugins/               # Plugins du CLI
│   ├── base/CLI_plugin.py     # Brique de base pour les plugins
│   ├── get.py, comp.py, ...   # Commandes
│
├── analysis/                  # Benchmarks/Calculs sur les données
│   ├── clusters.py            # Clustering hiérarchique
│   ├── comp_prediction.py     # Benchmark de prédiction (ESOL)
│   ├── comp_similarity.py     # Benchmark de continuité du voisinage(ESOL)
│   ├── purity_analysis.py     # Etant donné une famille ontologique, regarde les clusters
│   ├── dominante_families.py  # Evalue la constitution des clusters
│   └── correlation.py         # Corrélation syntaxique vs sémantique
│
├── similarites/               #
│   ├── cwl_kernel.py          # Cellular Weisfeiler-Lehman
│   ├── builtin_similarity.py  # Fingerprints RDKit (Morgan, etc.)
│   └── ontology_similarity.py # Similarité basée sur l'arbre ChEBI
│
├── Iso/                       # Isomorphisme de graphes
│
├── datasets/                  # Stockage des jeux de données externes
│   └── esol.csv               # Dataset ESOL (téléchargé automatiquement)
│
└── DB_updates/                # Métadonnées pour le maintien à jour des BDD
```

### Notes

- **`cli_plugins/`** : Le CLI détecte automatiquement les commandes, tant que le fichier contient une classe qui implemente CLIPlugin (qui se trouve dans cli_plugins/base). Chacunes de ces classes définit une commande CLI (ex: `get`, `comp`, `clustering`). Cela rend le projet très facile à étendre.
- **`similarites/`** : Chaque similarité doit implémenter soit la brique de base BaseGraphSimilarity (dans le cas de noyaux personnalisés), soit BaseSimilarityFromFingerprint (dans le cas de fingerprints).

---

## Installation

### 1) Pré-requis

- **Python 3.10+**
- Accès internet (pour télécharger la base ChEBI et le dataset ESOL)

### 3) Dépendances

Le projet s'appuie sur des bibliothèques scientifiques standards en Python :

```bash
pip install rdkit numpy scipy scikit-learn tqdm requests matplotlib
```

---

## Démarrage Rapide (Quickstart)

Le point d'entrée unique est le script `mol_cli.py`.

IMPORTANT : https://www.ebi.ac.uk/chebi/CHEBI:id ou CHEBI:id ou id seul ne font pas de différence dans les commandes du CLI.

**1. Préparer la base de données locale (Téléchargement initial)**

```bash
python3 ./mol_cli.py update
```

**2. Récupérer une molécule spécifique depuis ChEBI**

```bash
python3 ./mol_cli.py get --chebi-id 15377
```

**3. Comparer la similarité entre 2 molécules**

```bash
python3 ./mol_cli.py comp --id1 100 --id2 101 --fingerprint cwl --method tanimoto
```

**4. Lancer un benchmark complet sur le dataset ESOL**

```bash
python3 ./mol_cli.py comp_prediction --operation all --dataset datasets/esol.csv --fingerprints morgan,rdkit,cwl --models ridge,rf --similarities tanimoto,dice,cosine --folds 5
```

---

## Commandes CLI Détaillées

IMPORTANT : https://www.ebi.ac.uk/chebi/CHEBI:id ou CHEBI:id ou id seul ne font pas de différence dans les commandes du CLI.

### -h

Afin d'avoir la liste des commandes

```bash
python3 ./mol_cli.py -h
```

Pour chaque commande, pour avoir des informations sur ses arguments :

```bash
python3 ./mol_cli command -h
```

`get`

Récupère une molécule depuis l'API ChEBI et l'insère dans la base locale SQLite.

```bash
python3 ./mol_cli.py get --chebi-id 15377
```

### `update`

Vérifie si la base de données locale est à jour par rapport aux serveurs ChEBI et télécharge les mises à jour si nécessaire.

```bash
python3 ./mol_cli.py update
```

### `iso`

Teste l'isomorphisme exact entre deux molécules présentes dans la base locale.

```bash
python3 ./mol_cli.py iso --id1 100 --id2 101
```

### `comp`

Calcule le score de similarité entre deux molécules.

```bash
python3 ./mol_cli.py comp --id1 100 --id2 101 --fingerprint morgan --method tanimoto
```

*Options :*

- `--fingerprint` : `cwl`, `morgan`, `rdkit`
- `--method` : `tanimoto`, `dice`, `cosine`

### `clustering`

Effectue un clustering hiérarchique sur les molécules de la base et  options d'analyse.
*Opérations disponibles : `run`, `load`, `purity`, `dominant`, `cdf`.*

**Exemple : Lancer le clustering**

```bash
python3 ./mol_cli.py clustering --operation run --kernel cwl --similarity tanimoto --threshold 0.75 --sample-size 2000 --output-file clusters_data.json
```

**Exemple : Analyser la pureté d'une famille**

```bash
python3 ./mol_cli.py clustering --operation purity --input-file clusters_data.json --family flavonoids
```

### `comp_prediction`

Lance des benchmarks de prédiction de propriétés (logS) ou de similarité sur le dataset ESOL.
*Opérations disponibles : `prediction`, `similarity`, `all`.*

**Exemple : Prédiction seule**

```bash
python3 ./mol_cli.py comp_prediction --operation prediction --dataset datasets/esol.csv --fingerprints morgan,rdkit,cwl --models ridge,rf --folds 5
```

*(Note : Si le fichier `esol.csv` n'existe pas, il sera téléchargé automatiquement).*

### `correlation`

Lance un benchmark pour évaluer la corrélation entre la similarité syntaxique (structurelle) et sémantique (ontologique) sur un sous-ensemble de familles.

```bash
python3 ./mol_cli.py correlation --per-family 40 --kernel-similarity cosine --fingerprint-method cwl
```
