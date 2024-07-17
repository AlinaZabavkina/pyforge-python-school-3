from rdkit import Chem

smiles_container = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
smiles_structure = "c1ccccc1"

def substructure_search(smiles_container, smiles_structure):
    molecule_list = [Chem.MolFromSmiles(x) for x in smiles_container]
    substructure = Chem.MolFromSmiles(smiles_structure)
    matches = [Chem.MolToSmiles(x) for x in molecule_list if x.HasSubstructMatch(substructure)]
    return matches

print(substructure_search(smiles_container, smiles_structure))
