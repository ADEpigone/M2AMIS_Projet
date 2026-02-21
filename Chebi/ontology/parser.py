from Chebi.ontology.ontology_tree import OntologyTree, OntologyNode


def parse_obo(obo_text: str) -> OntologyTree:
    """
    Parse un fichier OBO (ChEBI) et construit un OntologyTree.
    
    Format OBO : https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html
    """
    tree = OntologyTree()
    
    current_term: dict = None
    in_term = False

    for line in obo_text.splitlines():
        line = line.strip()

        if line == "[Term]":
            # Sauvegarder le terme précédent
            if current_term and current_term.get("id"):
                _commit_term(tree, current_term)
            current_term = {}
            in_term = True
            continue

        if line == "[Typedef]" or line == "":
            if in_term and current_term and current_term.get("id"):
                _commit_term(tree, current_term)
                current_term = None
            if line == "[Typedef]":
                in_term = False
            continue

        if not in_term or current_term is None:
            continue

        if ": " not in line:
            continue

        key, _, value = line.partition(": ")
        key = key.strip()
        value = value.strip()

        if key == "id":
            current_term["id"] = value
        elif key == "name":
            current_term["name"] = value
        elif key == "def":
            # def: "blabla" [ref]
            current_term["definition"] = value.split('"')[1] if '"' in value else value
        elif key == "is_a":
            # is_a: CHEBI:12345 ! some name
            parent_id = value.split("!")[0].strip()
            current_term.setdefault("is_a", []).append(parent_id)
        elif key == "relationship":
            # relationship: has_role CHEBI:12345
            parts = value.split("!")  [0].strip().split(None, 1)
            if len(parts) == 2:
                rel_type, target = parts
                current_term.setdefault("relationships", []).append((rel_type, target))
        elif key == "synonym":
            # synonym: "blabla" EXACT [...]
            if '"' in value:
                syn = value.split('"')[1]
                current_term.setdefault("synonyms", []).append(syn)
        elif key == "is_obsolete" and value == "true":
            current_term = None
        elif key == "property_value":
            # property_value: http://purl.obolibrary.org/obo/chebi/formula "C2H6O" xsd:string
            current_term.setdefault("property_values", []).append(value)

    # Dernier terme
    if current_term and current_term.get("id"):
        _commit_term(tree, current_term)

    tree.build_roots()
    return tree


def _commit_term(tree: OntologyTree, term: dict):
    """Ajoute un terme parsé dans l'arbre."""
    chebi_id = term["id"]
    node = tree.get_or_create_node(chebi_id)
    node.name = term.get("name", "")
    node.definition = term.get("definition", "")

    # Propriétés OBO (formula, mass, etc.)
    for pv in term.get("property_values", []):
        _parse_property_value(node, pv)

    # is_a -> parents
    for parent_id in term.get("is_a", []):
        parent_node = tree.get_or_create_node(parent_id)
        node.add_parent(parent_node)

    # relationships
    for rel_type, target_id in term.get("relationships", []):
        node.add_relationship(rel_type, target_id)
        if rel_type in {"has_role", "RO:0000087"}:
            node.add_role(target_id)


def _parse_property_value(node: OntologyNode, pv_str: str):
    """
    Parse une ligne property_value OBO.
    Ex: http://purl.obolibrary.org/obo/chebi/formula "C2H6O" xsd:string
    
    Propriétés extraites : formula, mass, charge, monoisotopicMass, smiles, inchi, inchikey
    """
    parts = pv_str.split('"')
    if len(parts) >= 2:
        prop_uri = parts[0].strip()
        prop_value = parts[1]
        # Extraire le nom court de la propriété
        prop_name = prop_uri.rsplit("/", 1)[-1] if "/" in prop_uri else prop_uri
        
        # Convertir les valeurs numériques
        if prop_name in ['mass', 'monoisotopicMass', 'charge']:
            try:
                prop_value = float(prop_value)
            except (ValueError, TypeError):
                pass
        
        node.add_relationship(prop_name, prop_value)