import os
import sys
import time
import faiss
import pickle
import psycopg2
import numpy as np
import psycopg2.extras


batch_size = int(sys.argv[1])
output_dir = 'simulation_files/batch_indexes/'  
batch_state_path = 'simulation_files/batch_state.txt'
write_frequency = 10    # 1*batch_size will be size of one faiss file
# gpu_resources = faiss.StandardGpuResources()

class LSHIndex:
    def __init__(self, dimension, num_bits=1024, use_gpu=False):
        self.dimension = dimension
        self.num_bits = num_bits
        self.use_gpu = use_gpu
        self.create_index()

    def create_index(self):
        self.base_index = faiss.IndexLSH(self.dimension, self.num_bits)
        if self.use_gpu:
            self.res = gpu_resources
            self.base_index = faiss.index_cpu_to_gpu(self.res, 0, self.base_index)
        self.index = faiss.IndexIDMap(self.base_index)

    def add(self, vectors, labels):
        self.index.add_with_ids(np.array(vectors, dtype=np.float32), np.array(labels, dtype=np.int64))

    def write_index(self, path):
        faiss.write_index(self.index, path)
        self.create_index()  

def save_last_processed_index(last_index):
    with open(batch_state_path, 'w') as f:
        f.write(str(last_index))

def load_last_processed_index():
    if os.path.exists(batch_state_path):
        with open(batch_state_path, 'r') as f:
            content = f.read()
            if content.strip():
                return int(content)
    return 0

def process_rows(rows):
    vectors = []
    ids = []
    for row in rows:
        if row[1] is not None:
            try:
                vector = pickle.loads(row[1])
                vectors.append(vector)
                ids.append(row[0])
            except Exception as e:
                print(f"Error processing row {row[0]}: {e}")
    if vectors:
        vectors = np.vstack(vectors)
    else:
        vectors = np.empty((0, 1024))  
    ids = np.array(ids, dtype=np.int64)
    return vectors, ids

conn = psycopg2.connect(host="db", dbname="fingerprint_db", user="kailash", password="enamine", port=5432)
conn.autocommit = False  
cursor = conn.cursor(name='server_cursor')

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

last_processed_index = load_last_processed_index()
lsh_index = LSHIndex(dimension=1024, num_bits=1024, use_gpu=False)

print("Connecting to the database and initializing index...")
batch_count = last_processed_index // batch_size
print(f"Indexing will start from Batch {batch_count + 1}")

cursor.execute(f"SELECT db_index, fingerprints FROM fingerprints_table WHERE db_index > {last_processed_index} ORDER BY db_index ASC")

last_index_used = last_processed_index

try:
    while True:
        start_time = time.time()
        rows = cursor.fetchmany(batch_size)
        if not rows:
            break
        batch_vectors, batch_ids = process_rows(rows)
        lsh_index.add(batch_vectors, batch_ids)
        last_index_used = batch_ids[-1]  

        if (batch_count + 1) % write_frequency == 0:
            batch_index_path = os.path.join(output_dir, f'lsh_index_batch_{((batch_count // write_frequency) % 10) + 1}.faiss')
            lsh_index.write_index(batch_index_path)
            save_last_processed_index(last_index_used)
            print(f"Index file writing to: {batch_index_path}")
        batch_count += 1
        end_time = time.time()
        print(f"Batch {batch_count} processed in {(end_time - start_time) / 60:.2f} min")
finally:
    cursor.close()
    conn.close()

print("Index generation completed.")
