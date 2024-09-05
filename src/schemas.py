from pydantic import BaseModel


class MoleculeCreate(BaseModel):
    # molecule_id: int
    smiles_structure: str

class MoleculeUpdate(BaseModel):
    smiles_structure: str

