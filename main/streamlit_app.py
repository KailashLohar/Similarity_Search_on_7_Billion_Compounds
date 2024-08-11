import os
import time
import faiss
import psycopg2
import numpy as np
import pandas as pd
import streamlit as st
from rdkit import Chem
from concurrent.futures import ThreadPoolExecutor, as_completed
from rdkit.Chem import rdFingerprintGenerator, DataStructs

def query_index(file_path, dense_fp, k):
    lsh_index = faiss.read_index(file_path)
    _, query_indices = lsh_index.search(dense_fp.reshape(1, -1), k)
    return query_indices[0]

def compute_tanimoto_similarity(query_mol, mol):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    fp1 = gen.GetFingerprint(query_mol)
    fp2 = gen.GetFingerprint(mol)
    tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)
    return round(tanimoto, 4)

def process_single_index(file_path, dense_fp, k, query_mol, cursor, batch_num, total_batches):
    step_start_time = time.time()
    
    query_indices = query_index(file_path, dense_fp, k)
    if query_indices.size > 0:
        query_indices = [int(idx) for idx in query_indices]
        placeholders = ','.join(['%s'] * len(query_indices))
        cursor.execute(f'SELECT id, smiles FROM fingerprints_table WHERE db_index IN ({placeholders})', tuple(query_indices))
        results = cursor.fetchall()

        df = pd.DataFrame(results, columns=['id', 'smiles'])
        df['Tanimoto_similarity'] = df['smiles'].apply(lambda x: compute_tanimoto_similarity(query_mol, Chem.MolFromSmiles(x)))
        df = df.sort_values(by='Tanimoto_similarity', ascending=False).head(k).reset_index(drop=True)
        df = df.rename(columns={'id': 'ID', 'smiles': 'SMILES', 'Tanimoto_similarity': 'Score'})
        
        step_end_time = time.time()
        total_step_time = int(step_end_time - step_start_time)
        st.sidebar.write(f"Processed batch {batch_num}/{total_batches} in {total_step_time} s")
        
        return df
    return pd.DataFrame()

def process_similarity_query(query_smiles: str, k: int, conn):
    query_mol = Chem.MolFromSmiles(query_smiles)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    query_fp = gen.GetFingerprint(query_mol)

    dense_fp = np.zeros((1024,), dtype=np.float32)
    Chem.DataStructs.ConvertToNumpyArray(query_fp, dense_fp)

    index_files = [os.path.join('simulation_files/batch_indexes', f) for f in os.listdir('simulation_files/batch_indexes') if f.endswith('.faiss')]
    total_batches = len(index_files)
    all_results = pd.DataFrame()

    with conn.cursor() as cursor:
        for batch_num, file_path in enumerate(index_files, start=1):
            result = process_single_index(file_path, dense_fp, k, query_mol, cursor, batch_num, total_batches)
            if not result.empty:
                all_results = pd.concat([all_results, result], ignore_index=True)

    return all_results

st.markdown("<h1 style='color: #8b2721;'>Enamine Database Explorer</h1>", unsafe_allow_html=True)
st.markdown("<h2 style='color: #0b6098; font-size: 22px;'>Similarity Search Across 6.7 Billion Compounds</h2>", unsafe_allow_html=True)

st.sidebar.image("logo.svg", width=200)
st.sidebar.title("Compound Similarity Query")
query_smiles = st.sidebar.text_input("Enter SMILES of the Query Compound:", "Cc1cnc2c(cccc2c1)S(=O)(=O)N1CCN(CC1)C(=O)Nc1ccc(F)cc1")
k = st.sidebar.number_input("Number of Compounds to Retrieve:", value=20, min_value=1)
submit_button = st.sidebar.button("Find Similar Compounds")

if submit_button:
    with psycopg2.connect(host="db", dbname="fingerprint_db", user="kailash", password="enamine", port=5432) as conn:
        index_files = [os.path.join('simulation_files/batch_indexes', f) for f in os.listdir('simulation_files/batch_indexes') if f.endswith('.faiss')]
        total_batches = len(index_files)
        st.sidebar.write(f"Total {total_batches} files")

        overall_start_time = time.time()
        all_results = process_similarity_query(query_smiles, k, conn)

        all_results = all_results.sort_values(by='Score', ascending=False).head(k).reset_index(drop=True)
        
        overall_end_time = time.time()
        elapsed_time = overall_end_time - overall_start_time

        st.write(f"Top {k} Similar Compounds")
        st.write(all_results)
        st.write(f"Total Execution Time: {elapsed_time:.2f} seconds")