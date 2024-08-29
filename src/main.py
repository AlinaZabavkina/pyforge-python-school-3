from fastapi import FastAPI, HTTPException, Depends, status, UploadFile, Query
from sqlalchemy.orm import Session
from rdkit import Chem
import logging
import csv
import io
from models import SessionLocal, engine, Molecule
from schemas import MoleculeUpdate, MoleculeCreate

app = FastAPI()

# Set up logging configuration
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Dependency to get the SQLAlchemy session
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

@app.get("/molecules/all", tags=["Molecules"], summary="Get all molecules")
def get_molecules(limit: int = Query(None, description="Limit the number of molecules returned"), db: Session = Depends(get_db)):
    """Get List all molecules. Endpoint will return a limited number of molecules based on the provided limit"""
    logger.info(f"Fetching up to {limit} molecules")
    query = db.query(Molecule)
    if limit:
        query = query.limit(limit)

    molecules = query.all()
    if molecules:
        logger.info(f"Fetched {len(molecules)} molecules")
        return (molecule for molecule in molecules)
    else:
        logger.warning("No molecules found")
        raise HTTPException(status_code=404, detail="Molecules are not found")


@app.get("/molecules/{molecule_id}", tags=["Molecules"], summary="Find molecule by id")
def get_molecule(molecule_id: int, db: Session = Depends(get_db)):
    """Get molecule by identifier. Endpoint will return molecule with searched id"""
    logger.info(f"Fetching molecule with ID: {molecule_id}")
    molecule = db.query(Molecule).filter(Molecule.molecule_id == molecule_id).first()
    if molecule:
        logger.info(f"Molecule with ID {molecule_id} found")
        return molecule
    logger.warning(f"Molecule with ID {molecule_id} not found")
    raise HTTPException(status_code=404, detail="Molecule is not found")

@app.get("/molecules", tags=["Molecules"], summary="Find molecule by substructure")
def substructure_search(smiles_structure: str, db: Session = Depends(get_db)):
    """Substructure search for all added molecules. Endpoint will return all structures that contain the searched substructure"""
    logger.info(f"Substructure search for: {smiles_structure}")
    substructure = Chem.MolFromSmiles(smiles_structure)
    if substructure is None:
        logger.error("Invalid SMILES substructure provided")
        raise HTTPException(status_code=400, detail="Invalid SMILES substructure")

    molecules = db.query(Molecule).all()
    list_of_matches = []
    for molecule in molecules:
        mol_structure = Chem.MolFromSmiles(molecule.smiles_structure)
        if mol_structure is None:
            logger.error(f"Invalid SMILES structure found in database: {molecule.smiles_structure}")
            raise HTTPException(status_code=400, detail=f"Invalid SMILES structure: {molecule.smiles_structure}")
        elif mol_structure.HasSubstructMatch(substructure):
            list_of_matches.append({"molecule_id": molecule.molecule_id, "smiles_structure": molecule.smiles_structure})

    if list_of_matches:
        logger.info(f"Found {len(list_of_matches)} matches for the substructure")
        return list_of_matches
    else:
        logger.warning("No molecules found with the given substructure")
        raise HTTPException(status_code=404, detail="No molecules found with the given substructure")

@app.put("/molecules/{molecule_id}", tags=["Molecules"], summary="Update molecule structure", response_description="Molecule was updated successfully")
def update_molecule(molecule_id: int, updated_molecule: MoleculeUpdate, db: Session = Depends(get_db)):
    """Updating a molecule by identifier. Endpoint will return updated structures"""
    logger.info(f"Updating molecule with ID: {molecule_id}")
    molecule = db.query(Molecule).filter(Molecule.molecule_id == molecule_id).first()
    if molecule:
        molecule.smiles_structure = updated_molecule.smiles_structure
        db.commit()
        db.refresh(molecule)
        logger.info(f"Molecule with ID {molecule_id} updated successfully")
        return molecule
    logger.warning(f"Molecule with ID {molecule_id} not found for update")
    raise HTTPException(status_code=404, detail="Molecule with searched id is not found")

@app.post("/molecules", status_code=status.HTTP_201_CREATED, tags=["Molecules"], summary="Create new molecule", response_description="Molecule was created successfully")
def add_molecule(new_molecule: MoleculeCreate, db: Session = Depends(get_db)):
    """Add molecule and its identifier. Endpoint will return newly created structure"""
    logger.info("Adding a new molecule")
    molecule = Molecule(smiles_structure=new_molecule.smiles_structure)
    db.add(molecule)
    db.commit()
    db.refresh(molecule)
    logger.info(f"New molecule added with ID: {molecule.molecule_id}")
    return molecule

@app.post("/uploadfile", tags=["File upload"], summary="Upload file with molecules", response_description="File is uploaded successfully")
async def upload_file(file: UploadFile, db: Session = Depends(get_db)):
    """Upload file with molecules. Endpoint will parse data in the file into dict and return full list of structures"""
    logger.info(f"Uploading file: {file.filename}")
    if file.filename.endswith('.csv'):
        contents = await file.read()
        csv_file = io.StringIO(contents.decode('utf-8'))
        reader = csv.reader(csv_file)
        next(reader)  # Skip the header row
        for molecule_id, mol in reader:
            structure = Chem.MolFromSmiles(mol)
            if structure is None:
                logger.error(f"Invalid SMILES substructure in file: {mol}")
                raise HTTPException(status_code=400, detail=f"Invalid SMILES substructure in your file: {mol}")
            else:
                molecule_instance = Molecule(smiles_structure=mol)
                db.add(molecule_instance)
        db.commit()
        molecules = db.query(Molecule).all()
        logger.info("File processed and molecules added to the database")
        return molecules
    else:
        logger.warning("Invalid file extension")
        raise HTTPException(status_code=404, detail="File extension is not allowed. Please upload csv file.")

@app.delete("/molecules/{molecule_id}", tags=["Molecules"], summary="Delete molecule by id", response_description="Molecule was deleted successfully")
def delete_molecule(molecule_id: int, db: Session = Depends(get_db)):
    """Delete a molecule by identifier. Endpoint will return deleted structure"""
    logger.info(f"Deleting molecule with ID: {molecule_id}")
    molecule = db.query(Molecule).filter(Molecule.molecule_id == molecule_id).first()
    if molecule:
        db.delete(molecule)
        db.commit()
        logger.info(f"Molecule with ID {molecule_id} deleted successfully")
        return molecule
    logger.warning(f"Molecule with ID {molecule_id} not found for deletion")
    raise HTTPException(status_code=404, detail="Molecule is not found")
