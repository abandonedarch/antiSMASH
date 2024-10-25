
# antismash4.py
import os
import json
import csv
from Bio import SeqIO

# AminoacidCounter Class
class AminoacidCounter:
    def __init__(self):
        self.aminoacids = {}

    def merge_cluster_aminoacids(self, cluster_aminoacids):
        for aminoacid_name, aminoacid_database_bundle in cluster_aminoacids.items():
            if aminoacid_name not in self.aminoacids:
                self.aminoacids[aminoacid_name] = {}
            local_bundle = self.aminoacids[aminoacid_name]
            for db, count in aminoacid_database_bundle.items():
                if db not in local_bundle:
                    local_bundle[db] = 0
                local_bundle[db] += count

# BufferWriter Class
class BufferWriter:
    DEFAULT_CAPACITY = 16 * 1024

    def __init__(self, path):
        self.path = path
        self.capacity = BufferWriter.DEFAULT_CAPACITY
        self.data = ""

    def write(self, data):
        self.data += data
        if len(self.data) >= self.capacity:
            self.flush()

    def flush(self):
        with open(self.path, 'a', encoding='utf-8') as f:
            f.write(self.data)
        self.data = ""

# DbBundle and Cluster Classes
with open('./aminoacidNames.json', 'r', encoding='utf-8') as f:
    AMINOACID_NAMES = json.load(f)

class DbBundle:
    def __init__(self):
        self.data = {
            "Stachelhaus code": 0,
            "pHMM": 0,
            "NRPSPredictor3 SVM": 0,
            "SANDPUMA ensemble": 0,
        }

    def gather_information(self, info, rename_behavior=None):
        rename_behavior = rename_behavior or (lambda x: x)
        for key, value in self.data.items():
            info[rename_behavior(key)] = value

class Cluster:
    SMCOG_OFFSET = 1000
    SMCOG_LENGTH = 301

    def __init__(self, name, directory_name):
        self.name = name
        self.cluster_type = ""
        self.directory_name = directory_name
        self.l_aminoacids = {}
        self.d_aminoacids = {}
        self.epimerization = 0
        self.smcog_arr = [0] * Cluster.SMCOG_LENGTH
        self.starter_d_asn = DbBundle()
        self.starter_l_asn = DbBundle()
        self.serpro = DbBundle()

        Cluster.import_aminoacid_names(self.l_aminoacids, AMINOACID_NAMES['l'])
        Cluster.import_aminoacid_names(self.d_aminoacids, AMINOACID_NAMES['d'])

    @staticmethod
    def import_aminoacid_names(target_property, hardcoded_property):
        for name in hardcoded_property:
            target_property[name] = DbBundle()

    def gather_information(self):
        info = {
            'fileName': self.directory_name,
            'clusterName': self.name,
            'clusterType': self.cluster_type,
            'epimerization': self.epimerization,
        }
        Cluster.gather_information_for_aminoacids(info, "l", self.l_aminoacids)
        Cluster.gather_information_for_aminoacids(info, "d", self.d_aminoacids)
        self.starter_l_asn.gather_information(info, lambda db_name: f"{db_name}_starter_l_asn")
        self.starter_d_asn.gather_information(info, lambda db_name: f"{db_name}_starter_d_asn")
        self.serpro.gather_information(info, lambda db_name: f"{db_name}_ser_pro")
        for i in range(Cluster.SMCOG_LENGTH):
            info[f"SMCOG{i + Cluster.SMCOG_OFFSET}"] = self.smcog_arr[i]
        return info

    @staticmethod
    def gather_information_for_aminoacids(info, prefix, self_property):
        for aminoacid_name, aminoacid_bundle in self_property.items():
            aminoacid_bundle.gather_information(info, lambda db_name: f"{db_name}_{prefix}_{aminoacid_name}")

# Genus Class
class Genus:
    def __init__(self, name):
        self.name = name
        self.a_cluster = []

    def read(self, geneclusters, details_data, directory_name):
        for cluster_name, g_cluster in geneclusters.items():
            cluster = Cluster(cluster_name, directory_name)
            cluster.read(g_cluster, details_data.get(cluster_name, {}))
            self.a_cluster.append(cluster)

    def run_on_directory(self, path, directory_name):
        a_file_name = os.listdir(path)
        if 'geneclusters.js' in a_file_name:
            # Assuming geneclusters and details_data are defined in geneclusters.js
            geneclusters = {}
            details_data = {}
            self.read(geneclusters, details_data, directory_name)
            return
        for file_name in a_file_name:
            full_path = os.path.join(path, file_name)
            if os.path.isdir(full_path):
                self.run_on_directory(full_path, directory_name)

# Extract antiSMASH Structures Function for antiSMASH 4 using GBK files
def extract_antismash4_structures(gbk_file):
    clusters_data = []

    # Parse GBK file using Biopython
    with open(gbk_file, 'r') as file:
        for record in SeqIO.parse(file, 'genbank'):
            for feature in record.features:
                if feature.type == 'region' and 'product' in feature.qualifiers:
                    cluster_info = {
                        'cluster_number': feature.qualifiers.get('region_number', ['N/A'])[0],
                        'cluster_type': feature.qualifiers.get('product', ['N/A'])[0]
                    }

                    # Extract chemical structure information if available
                    if 'chemical_structure' in feature.qualifiers:
                        cluster_info['predicted_structure'] = feature.qualifiers['chemical_structure'][0]
                        cluster_info['smiles'] = feature.qualifiers.get('smiles', [''])[0]
                        cluster_info['molecular_formula'] = feature.qualifiers.get('molecular_formula', [''])[0]

                    clusters_data.append(cluster_info)

    return clusters_data

# Process antiSMASH 4 Results Function
def process_antismash4_results(directory):
    all_clusters_data = []

    # Loop through files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.gbk'):
            gbk_file = os.path.join(directory, filename)
            print(f"Processing: {gbk_file}")
            clusters_data = extract_antismash4_structures(gbk_file)
            all_clusters_data.extend(clusters_data)

    return all_clusters_data

# Save to CSV Function
def save_to_csv(data, output_file):
    keys = ['cluster_number', 'cluster_type', 'predicted_structure', 'smiles', 'molecular_formula']
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        writer.writeheader()
        writer.writerows(data)

    print(f"Data saved to {output_file}")

# Main Script
buffer_writer = BufferWriter('result_test.txt')
all_genus = {}

# Write empty JSON file
with open("./aminoacidNames.json", 'w', encoding='utf-8') as f:
    json.dump({'l': [], 'd': []}, f)

l_counter = AminoacidCounter()
d_counter = AminoacidCounter()

if os.name != 'nt':
    path = "/home/zzhong/Data/batch_output_vanco"
else:
    path = "E:/actinobacteria/Streptomyces"

def flush_genus():
    global all_genus
    for genus_name, genus in all_genus.items():
        for cluster in genus.a_cluster:
            info = cluster.gather_information()
            if not hasattr(flush_genus, 'has_title'):
                header_line = '	'.join(info.keys())
                buffer_writer.write(header_line + '
')
                flush_genus.has_title = True
            line = '	'.join(map(str, info.values()))
            buffer_writer.write(line + '
')
    all_genus = {}

# Process directories
for file_name in os.listdir(path):
    full_path = os.path.join(path, file_name)
    if os.path.isdir(full_path):
        genus_name = file_name.split('_')[0]
        genus = all_genus.get(genus_name, Genus(genus_name))
        genus.run_on_directory(full_path, file_name)
        all_genus[genus_name] = genus
        flush_genus()

buffer_writer.flush()

# Write results to JSON
with open("./aminoacidNames.json", 'w', encoding='utf-8') as f:
    json.dump({'l': list(l_counter.aminoacids.keys()), 'd': list(d_counter.aminoacids.keys())}, f)

# Process antiSMASH 4 results and save to CSV
antismash_results_dir = 'path_to_antismash_results_directory'  # Replace with your directory path
output_csv_file = 'extracted_structures.csv'
clusters_data = process_antismash4_results(antismash_results_dir)
save_to_csv(clusters_data, output_csv_file)
