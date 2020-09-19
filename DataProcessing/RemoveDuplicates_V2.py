# Changes over V1:
# 1. Removes duplicates from the Selected Biological Assemblies Folder

# Wacky error:
# Traceback (most recent call last):
#   File ".\RemoveDuplicates_V2.py", line 105, in <module>
#     if has_duplicate(path, pdb_list, max_sequence_similarity, max_structure_similarity):
#   File ".\RemoveDuplicates_V2.py", line 74, in has_duplicate
#     if is_equivalent(pdb, other_pdb, max_sequence_similarity, max_structure_similarity):
#   File ".\RemoveDuplicates_V2.py", line 32, in is_equivalent
#     structure_1 = parser_1.get_structure("structure_1", pdb1_filename)
#   File "C:\Users\alexp\Anaconda3\envs\Bioinformatics\lib\site-packages\Bio\PDB\MMCIFParser.py", line 64, in get_structure
#     self._build_structure(structure_id)
#   File "C:\Users\alexp\Anaconda3\envs\Bioinformatics\lib\site-packages\Bio\PDB\MMCIFParser.py", line 140, in _build_structure
#     chainid = chain_id_list[i]
# IndexError: list index out of range

from os.path import join
from os import listdir, remove
from Bio import SeqIO
import sys
import shutil
from Bio.PDB import MMCIFParser, PDBParser
import My_Library as ml
import numpy as np

from My_Library import *

def is_equivalent(pdb1_filename, pdb2_filename, max_sequence_similarity, max_structure_similarity):
    global cif_directory
    global sequence_cache
    global temp_path
    # Get the pdb_ids
    print(f"COMPARING {pdb1_filename} to {pdb2_filename}")
    pdb_id_1 = os.path.basename(pdb1_filename).split(".")[0]
    pdb_id_2 = os.path.basename(pdb2_filename).split(".")[0]
    # Set the parsers
    if pdb1_filename.split(".")[-1] == "cif":
        parser_1 = MMCIFParser(QUIET=True)
    elif pdb1_filename.split(".")[-1] in ("ent", "pdb"):
        parser_1 = PDBParser(QUIET=True)
    if pdb2_filename.split(".")[-1] == "cif":
        parser_2 = MMCIFParser(QUIET=True)
    elif pdb2_filename.split(".")[-1] in ("ent", "pdb"):
        parser_2 = PDBParser(QUIET=True)
    # Get the structures and get the sequences for each chain
    structure_1 = parser_1.get_structure("structure_1", pdb1_filename)
    sequences_1 = []
    for chain in structure_1.get_chains():
        chain_id = f"{pdb_id_1}{chain.get_id()}"
        if chain_id in sequence_cache: 
            # If we've found this chain before, get its sequence from cache
            sequences_1.append(sequence_cache[chain_id])
        else:
            # If not, get it, and add it to the cache
            seq = ml.load_chain_sequence(structure_1, chain.get_id())
            sequences_1.append(seq)
            sequence_cache[chain_id] = seq
    structure_2 = parser_2.get_structure("structure_2", pdb2_filename)
    sequences_2 = []
    for chain in structure_2.get_chains():
        chain_id = f"{pdb_id_2}{chain.get_id()}"
        if chain_id in sequence_cache: 
            # If we've found this chain before, get its sequence from cache
            sequences_2.append(sequence_cache[chain_id])
        else:
            # If not, get it, and add it to the cache
            seq = ml.load_chain_sequence(structure_2, chain.get_id())
            sequences_2.append(seq)
            sequence_cache[chain_id] = seq
    # Now, create a matrix that describes the PID of each combination of 
    # chains between structures 1 and 2
    equivalence_matrix = np.zeros((len(sequences_1), len(sequences_2)))
    for i in range(equivalence_matrix.shape[0]):
        for j in range(equivalence_matrix.shape[1]):
            print(f"Sequence 1: {sequences_1[i]}")
            print(f"Sequence 2: {sequences_2[j]}")
            equivalence_matrix[i, j] = ml.global_alignment_score(sequences_1[i], sequences_2[j]) == 1
    equivalence_matrix = np.matrix(equivalence_matrix)
    # If there is a 1-to-1 mapping of chains from structure 1 to structure 2, then the structures
    # are duplicated. Check this by seeing whether the matrix is bijective
    rank = np.linalg.matrix_rank(equivalence_matrix)
    assert equivalence_matrix.shape[0] == equivalence_matrix.shape[1]
    if  rank == equivalence_matrix.shape[0] and rank == equivalence_matrix.shape[1]:
        return True
    else: 
        return False

def has_duplicate(pdb, pdb_list, max_sequence_similarity, max_structure_similarity):
    for other_pdb in pdb_list:
        if is_equivalent(pdb, other_pdb, max_sequence_similarity, max_structure_similarity):
            return True
    return False

def get_filename(pdb_id):
    return f"pdb{pdb_id.lower()}.ent"
def get_pdbID(filename):
    return 

if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) != 4:
        print(f"\nUsage: $ python RemoveDuplicates.py <Species Root Directory> <Max Sequence Similarity> <Max Structure Identity>\n")
        sys.exit()
    root_directory = sys.argv[1]
    max_sequence_similarity = float(sys.argv[2])
    max_structure_similarity = float(sys.argv[3])

    # Remove duplicates
    print("1. REMOVING DUPLICATES")
    for species_folder in os.listdir(root_directory):
        print(f"{species_folder}")
        duplicates_removed_directory = os.path.join(root_directory, species_folder, "Duplicates_Removed")
        if not os.path.exists(duplicates_removed_directory): os.mkdir(duplicates_removed_directory)
        cif_directory = os.path.join(root_directory, species_folder, "Selected_Biological_Assemblies")
        temp_path = os.path.join(root_directory, species_folder, "_Temp")
        pdb_list = []
        sequence_cache = {}
        for filename in listdir(cif_directory):
            path = os.path.join(cif_directory, filename)
            if has_duplicate(path, pdb_list, max_sequence_similarity, max_structure_similarity):
                print(f"    Removed {filename}")
            else:
                pdb_list.append(path)
                destination = os.path.join(duplicates_removed_directory, filename)
                if not os.path.exists(destination):
                    shutil.copyfile(path, destination)
                print(f"    Saved {filename}")
            

        

    

    
