from os.path import join
from os import listdir, remove
from Bio import SeqIO
import sys
import shutil

from My_Library import *

def is_equivalent(pdb1_filename, pdb2_filename, max_sequence_similarity, max_structure_similarity):
    global cif_directory
    global sequence_cache
    global temp_path
    if pdb1_filename in sequence_cache: 
        sequence_1 = sequence_cache[pdb1_filename]
    else: 
        sequence_1 = load_sequence(pdb1_filename)
        sequence_cache[pdb1_filename] = sequence_1
    if pdb2_filename in sequence_cache: 
        sequence_2 = sequence_cache[pdb2_filename]
    else: 
        sequence_2 = load_sequence(pdb2_filename)
        sequence_cache[pdb2_filename] = sequence_2
    if global_alignment_score(sequence_1, sequence_2) >= max_sequence_similarity:
        # and calculate_TMScore(pdb1_filename, pdb2_filename, temp_path) >= max_structure_similarity
        return True
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
        cif_directory = os.path.join(root_directory, species_folder, "Clean_Structures")
        temp_path = os.path.join(root_directory, species_folder, "_Temp")
        pdb_list = []
        sequence_cache = {}
        for filename in listdir(cif_directory):
            path = os.path.join(cif_directory, filename)
            if has_duplicate(path, pdb_list, max_sequence_similarity, max_structure_similarity):
                print(f"    Removed {filename}")
            else:
                pdb_list.append(path)
        duplicates_removed_directory = os.path.join(root_directory, species_folder, "Duplicates_Removed")
        if not os.path.exists(duplicates_removed_directory): os.mkdir(duplicates_removed_directory)
        print('    Copying non-duplicates to Duplicates_Removed')
        for non_duplicate_path in pdb_list:
            non_duplicate_filename = os.path.basename(non_duplicate_path)
            destination = os.path.join(duplicates_removed_directory, non_duplicate_filename)
            if not os.path.exists(destination):
                shutil.copyfile(non_duplicate_path, destination)

        

    

    
