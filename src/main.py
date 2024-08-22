from fastapi import FastAPI, HTTPException, Depends, status, UploadFile
from sqlalchemy.orm import Session
from rdkit import Chem

from models import SessionLocal, engine, Molecule
from schemas import MoleculeUpdate, MoleculeCreate
import csv
import io

app = FastAPI()

# Dependency to get the SQLAlchemy session
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


# smiles_container = [{"molecule_id": 1, "smiles_structure": "c1cc(C)ccc1"},
#                     {"molecule_id": 2, "smiles_structure": "CCO"},
#                     {"molecule_id": 3, "smiles_structure": "CC(=O)O"},
#                     {"molecule_id": 4, "smiles_structure": "CC(=O)Oc1ccccc1C(=O)O"}]

@app.get("/molecules/all", tags=["Molecules"], summary="Get all molecules")
def get_molecules(db: Session = Depends(get_db)):
    """Get List all molecules. Endpoint will return all existed molecules"""
    molecules = db.query(Molecule).all()
    if molecules:
        return molecules
    else:
        raise HTTPException(status_code=404, detail="Molecules are not found")

@app.get("/molecules/{molecule_id}", tags=["Molecules"], summary="Find molecule by id")
def get_molecule(molecule_id: int, db: Session = Depends(get_db)):
    """Get molecule by identifier. Endpoint will return molecule with searched id"""
    # Convert list to dictionary for O(1) lookups
    molecule = db.query(Molecule).filter(Molecule.molecule_id == molecule_id).first()
    if molecule:
        return molecule
    raise HTTPException(status_code=404, detail="Molecule is not found")

@app.get("/molecules", tags=["Molecules"], summary="Find molecule by substructure")
def substructure_search(smiles_structure: str, db: Session = Depends(get_db)):
    """Substructure search for all added molecules. Endpoint will return all structures that contain the searched substructure"""
    substructure = Chem.MolFromSmiles(smiles_structure)
    if substructure is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES substructure")

    molecules = db.query(Molecule).all()
    list_of_matches = []
    for molecule in molecules:
        mol_structure = Chem.MolFromSmiles(molecule.smiles_structure)
        if mol_structure is None:
            raise HTTPException(status_code=400, detail=f"Invalid SMILES structure: {molecule.smiles_structure}")
        elif mol_structure.HasSubstructMatch(substructure):
            list_of_matches.append({"molecule_id": molecule.molecule_id, "smiles_structure": molecule.smiles_structure})

    if list_of_matches:
        return list_of_matches
    else:
        raise HTTPException(status_code=404, detail="No molecules found with the given substructure")

@app.put("/molecules/{molecule_id}", tags=["Molecules"], summary="Update molecule structure", response_description="Molecule was updated successfully")
def update_molecule(molecule_id: int, updated_molecule_id: MoleculeUpdate, db: Session = Depends(get_db)):
    """Updating a molecule by identifier.Endpoint will return updated structures"""
    molecule = db.query(Molecule).filter(Molecule.molecule_id == molecule_id).first()
    if molecule:
        molecule.smiles_structure = updated_molecule_id.smiles_structure
        db.commit()
        db.refresh(molecule)
        return molecule
    raise HTTPException(status_code=404, detail="Molecule with searched id is not found")

@app.post("/molecules", status_code =status.HTTP_201_CREATED, tags=["Molecules"], summary="Create new molecule",response_description="Molecule was created successfully")
def add_molecule(new_molecule: MoleculeCreate, db: Session = Depends(get_db)):
    """Add molecule and its identifier. Endpoint will return newly created structure"""
    molecule = Molecule(smiles_structure=new_molecule.smiles_structure)
    db.add(molecule)
    db.commit()
    db.refresh(molecule)
    return molecule

@app.post("/uploadfile", tags=["File upload"], summary="Upload file with molecules", response_description="File is uploaded successfully")
async def upload_file(file: UploadFile, db: Session = Depends(get_db)):
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
                raise HTTPException(status_code=400, detail=f"Invalid SMILES substructure in your file: {mol}")
            else:
                molecule_instance = Molecule(smiles_structure=mol)
                db.add(molecule_instance)
        db.commit()
        molecules = db.query(Molecule).all()
        return molecules
    else:
        raise HTTPException(status_code=404, detail="File extension is not allowed. Please upload csv file.")

@app.delete("/molecules/{molecules_id}", tags=["Molecules"], summary="Delete molecule by id",response_description="Molecule was deleted successfully")
def delete_molecule(molecule_id: int, db: Session = Depends(get_db)):
    """Delete a molecule by identifier. Endpoint will return deleted structure"""
    # Convert list to dictionary for O(1) lookups
    molecule = db.query(Molecule).filter(Molecule.molecule_id == molecule_id).first()
    if molecule:
        db.delete(molecule)
        db.commit()
        return molecule
    raise HTTPException(status_code=404, detail="Molecule is not found")

