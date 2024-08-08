from rdkit import Chem
from fastapi import FastAPI, status, HTTPException, UploadFile, Query, Body
from typing import List, Optional
import csv
import io

from src.models import MoleculeUpdate, MoleculeCreate

app = FastAPI()

smiles_container = [{"molecule_id": 1, "smiles_structure": "c1cc(C)ccc1"},
                    {"molecule_id": 2, "smiles_structure": "CCO"},
                    {"molecule_id": 3, "smiles_structure": "CC(=O)O"},
                    {"molecule_id": 4, "smiles_structure": "CC(=O)Oc1ccccc1C(=O)O"}]

@app.get("/molecules/all", tags=["Molecules"], summary="Get all molecules")
def get_molecules():
    """Get List all molecules. Endpoint will return all existed molecules"""
    if smiles_container:
        return smiles_container
    else:
        raise HTTPException(status_code=404, detail="Molecules are not found")

@app.get("/molecules/{molecule_id}", tags=["Molecules"], summary="Find molecule by id")
def get_molecule(molecule_id: int):
    """Get molecule by identifier. Endpoint will return molecule with searched id"""
    # Convert list to dictionary for O(1) lookups
    smiles_dict = {molecule["molecule_id"]: molecule for molecule in smiles_container}
    molecule = smiles_dict.get(molecule_id)
    if molecule:
        return molecule
    raise HTTPException(status_code=404, detail="Molecule is not found")

@app.get("/molecules", tags=["Molecules"], summary="Find molecule by substructure")
def substructure_search(smiles_structure: str, list_of_structures: Optional[List[MoleculeCreate]] = Body(None)):
    """Substructure search for all added molecules. Endpoint will return all structures that contain the searched substructure"""
    substructure = Chem.MolFromSmiles(smiles_structure)
    if substructure is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES substructure")
    list_of_matches = []
    for molecule in smiles_container:
        molecule_copy = molecule.copy()
        molecule_copy["smiles_structure"] = Chem.MolFromSmiles(molecule_copy["smiles_structure"])
        if molecule_copy["smiles_structure"] is None:
            raise HTTPException(status_code=400, detail=f"Invalid SMILES structure: {molecule}")
        elif molecule_copy["smiles_structure"].HasSubstructMatch(substructure):
            molecule_copy["smiles_structure"] = Chem.MolToSmiles(molecule_copy["smiles_structure"])
            list_of_matches.append(molecule_copy)

    if list_of_matches:
        return list_of_matches
    else:
        raise HTTPException(status_code=404, detail="No molecules found with the given substructure")

@app.put("/molecules/{molecule_id}", tags=["Molecules"], summary="Update molecule structure", response_description="Molecule was updated successfully")
def update_molecule(molecule_id: int, updated_molecule_id: MoleculeUpdate):
    """Updating a molecule by identifier.Endpoint will return updated structures"""
    for molecule in smiles_container:
        if molecule['molecule_id'] == molecule_id:
            molecule['smiles_structure'] = updated_molecule_id.smiles_structure
            return molecule
    raise HTTPException(status_code=404, detail="Molecule with searched id is not found")

@app.post("/molecules", status_code =status.HTTP_201_CREATED, tags=["Molecules"], summary="Create new molecule",response_description="Molecule was created successfully")
def add_molecule(new_molecule: MoleculeCreate):
    """Add molecule and its identifier. Endpoint will return newly created structure"""
    for molecule in smiles_container:  # check whether id is already existed
        if molecule["molecule_id"] == new_molecule.molecule_id:
            raise HTTPException(status_code=400, detail="Molecule ID already exists")

    smiles_container.append(new_molecule.dict())
    return new_molecule.dict()

@app.post("/uploadfile", tags=["File upload"], summary="Upload file with molecules", response_description="File is uploaded successfully")
async def upload_file(file: UploadFile):
    """Upload file with molecules.Endpoint will parse data in the file into dict and return full list of structures"""
    if file.filename.endswith('.csv'):
        contents = await file.read()
        csv_file = io.StringIO(contents.decode('utf-8'))
        reader = csv.reader(csv_file)
        next(reader)  # Skip the header row
        # for row in reader:
        for molecule_id, mol in reader:
            # Create an instance of the Molecule class
            structure = Chem.MolFromSmiles(mol)
            if structure is None:
                raise HTTPException(status_code=400, detail=f"Invalid SMILES substructure in your file: {row}")
            else:
                molecule_instance = MoleculeCreate(molecule_id=int(molecule_id), smiles_structure=mol)
                # Convert the instance to a dictionary
                molecule_dict = molecule_instance.dict()
                smiles_container.append(molecule_dict)
        return smiles_container
    else:
        raise HTTPException(status_code=404, detail="File extension is not allowed. Please upload csv file.")

@app.delete("/molecules/{molecules_id}", tags=["Molecules"], summary="Delete molecule by id",response_description="Molecule was deleted successfully")
def delete_molecule(molecule_id: int):
    """Delete a molecule by identifier. Endpoint will return deleted structure"""
    # Convert list to dictionary for O(1) lookups
    smiles_dict = {molecule["molecule_id"]: molecule for molecule in smiles_container}
    molecule = smiles_dict.get(molecule_id)
    if molecule:
        index = smiles_container.index(molecule)
        deleted_molecule = smiles_container.pop(index)
        return deleted_molecule
    raise HTTPException(status_code=404, detail="Molecule is not found")

