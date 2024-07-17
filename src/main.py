from rdkit import Chem

SMILES_container = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
SMILES_structure = "c1ccccc1"

def substructure_search(SMILES_container, SMILES_structure):
    molecule_list = [Chem.MolFromSmiles(x) for x in SMILES_container]
    substructure = Chem.MolFromSmiles(SMILES_structure)
    matches = [Chem.MolToSmiles(x) for x in molecule_list if x.HasSubstructMatch(substructure)]
    return matches

print(substructure_search(SMILES_container, SMILES_structure))
