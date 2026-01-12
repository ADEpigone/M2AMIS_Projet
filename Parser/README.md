Pour retrieve sur chebi il y a une nouvelle API :

- https://www.ebi.ac.uk/chebi/backend/api/docs/#/public/public_compound_retrieve

  - On y retrouve par exemple dans la réponse :
  - ```json
    "default_structure": {
      "id": 27932,
      "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
      "standard_inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
      "standard_inchi_key": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
      "wurcs": null,
      "is_r_group": false
    },
    ```

Pour parse le SMILES rdkit est utilisable par exemple : https://www.rdkit.org/docs/GettingStartedInPython.html

Puis pour convertir en graphe c'est assez direct.
