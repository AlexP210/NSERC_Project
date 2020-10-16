import os
import sys
import My_Library as ml
import shutil
import pandas as pd

def is_equivalent(pdb1_filename, pdb2_filename, max_sequence_similarity, max_structure_similarity):
    global redundant_HOMs_directory
    global sequence_cache
    global temp_directory
    global similar_chains_frame
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
    name_1 = os.path.basename(pdb1_filename).split(".")[0]
    name_2 = os.path.basename(pdb2_filename).split(".")[0]
    percent_identity = ml.global_alignment_score(sequence_1, sequence_2)
    if percent_identity >= max_sequence_similarity:
        if percent_identity < 1:
            nw_score, percent_identity = ml.calculate_percent_identity(name_1, pdb1_filename, name_2, pdb2_filename, os.path.join(temp_directory, "Global_Alignment.fasta"))
            tm = ml.calculate_TMScore(pdb1_filename, pdb2_filename, os.path.join(temp_directory, "TM_Align.txt"), alignment=os.path.join(temp_directory, "Global_Alignment.fasta")) >= max_structure_similarity
            similar_chains_frame = similar_chains_frame.append({
                "PDB_ID_1":name_1,
                "PDB_ID_2":name_2,
                "Chain_1_Length":len(sequence_1),
                "Chain_2_Length":len(sequence_2),
                "Length_Difference":abs(len(sequence_1) - len(sequence_2)),
                "Alignment_Score":nw_score,
                "PID":percent_identity,
                "TM":tm}, ignore_index=True)
        return True, similar_chains_frame
    return False, similar_chains_frame

def has_duplicate(pdb, pdb_list, max_sequence_similarity, max_structure_similarity):
    global similar_chains_frame
    for other_pdb in pdb_list:
        eq, similar_chains_frame = is_equivalent(pdb, other_pdb, max_sequence_similarity, max_structure_similarity)
        if eq:
            return True, similar_chains_frame
    return False, similar_chains_frame


if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) < 4:
        print(f"Usage: $ python {os.path.basename(__file__)} <Species Root Directory> <Max Sequence Similarity> <Max Structure Similarity>")
        sys.exit()
    root_directory = sys.argv[1]
    max_sequence_similarity = float(sys.argv[2])
    max_structure_similarity = float(sys.argv[3])

    sequence_cache = {}
    pdb_list = []
    nonredundant_HOMs_directory = os.path.join(root_directory, "..", "NonRedundant_HOMs")
    similar_chains_path = os.path.join(root_directory, "..", "SimilarHOMsComparisons.csv")
    similar_chains_frame = pd.DataFrame(columns=["PDB_ID_1","PDB_ID_2","PDB_1_Length","PDB_2_Length","Length_Difference","Alignment_Score","PID","TM"])
    if not os.path.exists(nonredundant_HOMs_directory): os.mkdir(nonredundant_HOMs_directory)
    # Go through the extracted chains for each species, and remove chains that are too similar to be random
    for species_folder in os.listdir(root_directory):
        print(species_folder)
        species_directory = os.path.join(root_directory, species_folder)
        temp_directory = os.path.join(species_directory, "_Temp")
        redundant_HOMs_directory = os.path.join(species_directory, "Clean_Structures")

        for cif_filename in os.listdir(redundant_HOMs_directory):
            print("     ", cif_filename)
            path = os.path.join(redundant_HOMs_directory, cif_filename)
            dup, similar_chains_frame = has_duplicate(path, pdb_list, max_sequence_similarity, max_structure_similarity)
            if dup:
                print("         Duplicated")
                pass
            else:
                print("         Original")
                pdb_list.append(path)
                shutil.copy(path, os.path.join(nonredundant_HOMs_directory, cif_filename))
    similar_chains_frame.to_csv(similar_chains_path)