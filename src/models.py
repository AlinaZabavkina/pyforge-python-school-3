from sqlalchemy import create_engine, Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
# from dotenv import load_dotenv
import os

# load_dotenv()
# SQLAlchemy setup
DATABASE_URL = "sqlite:///./test.db"

engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()

# SQLAlchemy ORM models
class Molecule(Base):
    __tablename__ = "molecules"

    molecule_id = Column(Integer, primary_key=True, index=True, autoincrement=True)
    smiles_structure = Column(String, index=True)

# Create tables
Base.metadata.create_all(bind=engine)
