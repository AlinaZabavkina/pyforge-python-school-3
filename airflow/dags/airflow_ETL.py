from datetime import timedelta

import boto3
import pandas as pd
from airflow import DAG
from airflow.operators.python import PythonOperator
from rdkit import Chem
from rdkit.Chem import Descriptors
from sqlalchemy.orm import Session
import os

from utils.dates import days_ago
from src.models import SessionLocal, Molecule

TABLE_NAME = "molecules"  # replace with your table name
S3_BUCKET_NAME = "hw-bucket-alinazzz777"
S3_FILE_PATH = "processed_data.xlsx"
access_key = os.getenv("AWS_ACCESS_KEY_ID")
secret_key = os.getenv("AWS_SECRET_ACCESS_KEY")


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

# 1. Extract Task
def extract_data(**kwargs):
    # Simulate extracting SMILES and related columns from the database
    # Replace this with actual DB connection and query logic
    db: Session = next(get_db())

    molecules = db.query(Molecule).all()
    df = pd.DataFrame(molecules)
    return df.to_dict()

# 2. Transform Task
def transform_data(**kwargs):
    ti = kwargs['ti']
    df_dict = ti.xcom_pull(task_ids='extract_data')
    df = pd.DataFrame(df_dict)

    # Initialize lists to hold the computed properties
    molecular_weights = []
    logP_values = []
    tpsa_values = []
    h_donors = []
    h_acceptors = []
    lipinski_pass = []

    molecules_dict = df[0]
    # Calculate molecular properties using RDKit
    for key, molecule in molecules_dict.items():  # Assuming df[0] contains Molecule objects
        smiles = molecule.smiles_structure  # Accessing the SMILES structure
        mol = Chem.MolFromSmiles(smiles)

        if mol is not None:
            molecular_weights.append(Descriptors.MolWt(mol))
            logP_values.append(Descriptors.MolLogP(mol))
            tpsa_values.append(Descriptors.TPSA(mol))
            h_donors.append(Descriptors.NumHDonors(mol))
            h_acceptors.append(Descriptors.NumHAcceptors(mol))

            # Apply Lipinski's rule of 5
            lipinski_pass.append(
                all([
                    Descriptors.MolWt(mol) < 500,
                    Descriptors.MolLogP(mol) < 5,
                    Descriptors.NumHDonors(mol) <= 5,
                    Descriptors.NumHAcceptors(mol) <= 10
                ])
            )
        else:
            molecular_weights.append(None)
            logP_values.append(None)
            tpsa_values.append(None)
            h_donors.append(None)
            h_acceptors.append(None)
            lipinski_pass.append(False)

    # Add the calculated properties to the DataFrame
    df['Molecular_Weight'] = molecular_weights
    df['logP'] = logP_values
    df['TPSA'] = tpsa_values
    df['H_Donors'] = h_donors
    df['H_Acceptors'] = h_acceptors
    df['Lipinski_Pass'] = lipinski_pass

    # Push transformed data to XCom for the next task
    kwargs['ti'].xcom_push(key='transformed_data', value=df.to_dict())
    return df.to_dict()

# 3. Save Task
def save_to_s3(**kwargs):
    ti = kwargs['ti']
    df_dict = ti.xcom_pull(task_ids='transform_data')
    df = pd.DataFrame(df_dict)

    # Save dataframe as .xlsx locally
    if not os.path.exists("/tmp/"):
        os.makedirs("/tmp/")

    output_file = "/tmp/processed_data.xlsx"
    df.to_excel(output_file, index=False)

    s3_client = boto3.client(
        's3',
        # access_key=access_key,
        # secret_key=secret_key,
        region_name='eu-north-1'  # Optional
    )
    try:
        # Upload the file to S3
        s3_client.upload_file(output_file, S3_BUCKET_NAME, S3_FILE_PATH)
        print(f"File  uploaded to successfully.")
    except Exception as e:
        print(f"Error uploading file: {e}")

    # Cleanup
    os.remove(output_file)

# Define DAG
default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'start_date': days_ago(1),
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 1,
    'retry_delay': timedelta(minutes=5),
}

dag = DAG(
    'smiles_processing_dag',
    default_args=default_args,
    description='A DAG to process SMILES data and save results to S3',
    schedule_interval=timedelta(days=1),  # Runs daily
)

# Define Tasks
extract_task = PythonOperator(
    task_id='extract_data',
    python_callable=extract_data,
    dag=dag,
)

transform_task = PythonOperator(
    task_id='transform_data',
    python_callable=transform_data,
    provide_context=True,
    dag=dag,
)

save_task = PythonOperator(
    task_id='save_to_s3',
    python_callable=save_to_s3,
    provide_context=True,
    dag=dag,
)

extract_task >> transform_task >> save_task