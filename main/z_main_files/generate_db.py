import os
import sys
import time
import gzip
import pickle
import psycopg2

from pyspark.sql import SparkSession
from pyspark.sql.functions import udf
from pyspark.sql.types import BinaryType, LongType
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

input_file_path = sys.argv[1]
batch_size = int(sys.argv[2])

total_chunks = 10
already_processed_index = (1-1)*batch_size    # put number x from file "simulation_files/x_temp.gz" in place of 1

spark = SparkSession.builder \
                    .appName("FingerprintDB") \
                    .config("spark.master", "local[15]") \
                    .config("spark.driver.memory", "20g") \
                    .config("spark.jars", "/spark/jars/postgresql-42.2.27.jar") \
                    .config("spark.network.timeout", "600s") \
                    .config("spark.executor.heartbeatInterval", "30s") \
                    .getOrCreate()

def mol_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
        mol_fp = gen.GetFingerprint(mol)
        return pickle.dumps(mol_fp)
    else:
        return None

generate_fingerprint_udf = udf(mol_fingerprint, BinaryType())
spark.udf.register("generate_fingerprint", mol_fingerprint, BinaryType())

conn = psycopg2.connect(host="db", dbname="fingerprint_db", user="kailash", password="enamine", port=5432)
cursor = conn.cursor()

cursor.execute("""DROP TABLE IF EXISTS fingerprints_table;
                  CREATE TABLE fingerprints_table (db_index BIGINT PRIMARY KEY,
                                                   id TEXT,
                                                   smiles TEXT,
                                                   fingerprints BYTEA);""")
conn.commit()

output_dir = 'simulation_files'
os.makedirs(output_dir, exist_ok=True)

with gzip.open(input_file_path, 'rt') as input_file:
    header = input_file.readline()  
    for _ in range(already_processed_index):
        input_file.readline()

    already_processed_chunks = already_processed_index // batch_size
    for i in range(already_processed_chunks + 1 , total_chunks + 1):
        start_time = time.time() 
        output_file_path = f'{output_dir}/{i}_temp.gz'
        with gzip.open(output_file_path, 'wt') as output_file:
            output_file.write(header)
            for _ in range(batch_size):
                line = input_file.readline()
                if not line:
                    break
                output_file.write(line)
        
        df = spark.read.option("header", "true").csv(f"simulation_files/{i}_temp.gz")
        df = df.withColumn("db_index", df["db_index"].cast(LongType()))
        df = df.withColumn("fingerprints", generate_fingerprint_udf("smiles")).repartition(20)
        df = df.orderBy("db_index")
        df.write \
          .format("jdbc") \
          .option("url", "jdbc:postgresql://db:5432/fingerprint_db") \
          .option("dbtable", "fingerprints_table") \
          .option("user", "kailash") \
          .option("password", "enamine") \
          .option("batchsize", batch_size) \
          .option("driver", "org.postgresql.Driver") \
          .mode("append") \
          .option("conflict", "append") \
          .save()

        os.remove(output_file_path)
        end_time = time.time()
        print(f"Chunk {i}/{total_chunks} processed in {(end_time - start_time)/60:.2f} min")

cursor.close()
conn.close()
spark.stop()