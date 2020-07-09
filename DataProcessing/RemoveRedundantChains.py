import os
import sys
import My_Library as ml
import shutil

def is_equivalent(pdb1_filename, pdb2_filename, max_sequence_similarity, max_structure_similarity):
    global redundant_chains_directory
    global sequence_cache
    global temp_directory
    if pdb1_filename in sequence_cache: 
        sequence_1 = sequence_cache[pdb1_filename]
    else: 
        sequence_1 = ml.load_sequence(pdb1_filename)
        sequence_cache[pdb1_filename] = sequence_1
    if pdb2_filename in sequence_cache: 
        sequence_2 = sequence_cache[pdb2_filename]
    else: 
        sequence_2 = ml.load_sequence(pdb2_filename)
        sequence_cache[pdb2_filename] = sequence_2
    if ml.global_alignment_score(sequence_1, sequence_2) >= max_sequence_similarity:
        # and calculate_TMScore(pdb1_filename, pdb2_filename, temp_path) >= max_structure_similarity
        return True
    return False

def has_duplicate(pdb, pdb_list, max_sequence_similarity, max_structure_similarity):
    for other_pdb in pdb_list:
        if is_equivalent(pdb, other_pdb, max_sequence_similarity, max_structure_similarity):
            return True
    return False


if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) < 4:
        print(f"Usage: $ python {os.path.basename(__file__)} <Species Root Directory> <Max Sequence Similarity> <Max Structure Similarity>")
        sys.exit()
    root_directory = sys.argv[1]
    max_sequence_similarity = float(sys.argv[2])
    max_structure_similarity = float(sys.argv[3])

    # Go through the extracted chains for each species, and remove chains that are too similar to be random
    for species_folder in os.listdir(root_directory):
        species_directory = os.path.join(root_directory, species_folder)
        redundant_chains_directory = os.path.join(species_directory, "Extracted_Chains")
        nonredundant_chains_directory = os.path.join(species_directory, "NonRedundant_Chains")
        if not os.path.exists(nonredundant_chains_directory): os.mkdir(nonredundant_chains_directory)

        sequence_cache = {}
        pdb_list = []
        for cif_filename in os.listdir(redundant_chains_directory):
            path = os.path.join(redundant_chains_directory, cif_filename)
            if has_duplicate(path, pdb_list, max_sequence_similarity, max_structure_similarity):
                pass
            else:
                pdb_list.append(path)
                shutil.copy(path, os.path.join(nonredundant_chains_directory, cif_filename))