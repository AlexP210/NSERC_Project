# Changes over V1:
#   1. calculates the rotation matrix
#   2. Excludes pairings corresponding to real heteromers
#   3. Does Ananas

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
import os.path
import os
import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import align, format_alignment
from My_Library import *
import random as rand
import pandas as pd
import sklearn as sk
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import My_Library as ml
import matplotlib.pyplot as plt
import numpy as np
import sklearn.decomposition as sk
from mpl_toolkits.mplot3d import Axes3D
import math as m

# def load_sequence(pdb_id, chain, directory):

#     # # for record in SeqIO.parse( os.path.join(directory, f"{pdb_id}.cif"), "cif-seqres" ):
#     # #     print(record.name)
#     # #     if record.id.split(":")[-1] == chain:
#     # #         return record.seq

#     parser = MMCIFParser(QUIET=True)
#     structure = parser.get_structure(f"{pdb_id}_{chain}", os.path.join(directory, f"{pdb_id}.cif"))
#     sequence = ""
#     for chain_structure in structure.get_chains():
#         if chain_structure.get_id() == chain:
#             for residue in chain_structure.get_residues():
#                 sequence += one_letter_code(residue.get_resname())
#             return sequence


if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) != 3:
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <Species Root Directory> <Number of Pairings>")
        sys.exit()
    root_directory = sys.argv[1]
    n_pairings = int(sys.argv[2])

    # Go over the list of chains, and compare any two chains that:
    #   1. Have both chains existing
    #   2. Were not previously compared
    #   3. Are not identical
    #   4. Have index as multiple of 500

    HOMs_directory = os.path.join(root_directory, "..", "NonRedundant_HOMs")
    output_file_path = os.path.join(root_directory, "..", "RandomComparisons.csv")
    random_comparisons = pd.DataFrame(columns=["PDB_ID_1","PDB_ID_2","PDB_1_Length","PDB_2_Length","Length_Difference","Alignment_Score","Alignment_Score_Adjusted","PID","TM"])

    temp_path = os.path.join(root_directory, "..", "_Temp")
    if not os.path.exists(temp_path): os.mkdir(temp_path)

    # Get the list of chains that we need to compare
    HOM_files = os.listdir(HOMs_directory)
    HOM_names = [cn.split(".")[0] for cn in HOM_files]

    # Prepare the output file

    compared = set()
    # Go through each combination of the chains
    for i in range(n_pairings):
        HOM_idx_a = rand.randint(0, len(HOM_names)-1)
        HOM_idx_b = rand.randint(0, len(HOM_names)-1)

        pdb_id_a = HOM_names[HOM_idx_a][:4].lower()
        pdb_a_filename = f"{pdb_id_a}.pdb"
        pdb_id_b = HOM_names[HOM_idx_b][:4].lower()
        pdb_b_filename = f"{pdb_id_b}.pdb"

        # If we have a bad choice, pick a new one
        while HOM_idx_a == HOM_idx_b or f"{pdb_id_a}{pdb_id_b}" in compared or f"{pdb_id_b}{pdb_id_a}" in compared:
            HOM_idx_a = rand.randint(0, len(HOM_names)-1)
            HOM_idx_b = rand.randint(0, len(HOM_names)-1)

            pdb_id_a = HOM_names[HOM_idx_a][:4].lower()
            pdb_a_filename = f"{pdb_id_a}.pdb"
            pdb_id_b = HOM_names[HOM_idx_b][:4].lower()
            pdb_b_filename = f"{pdb_id_b}.pdb"        
        compared.add(f"{pdb_id_a}{pdb_id_b}")

        print(f"Processing: {pdb_id_a} vs {pdb_id_b} # {i + 1}")
        # Get the paths to the chains, and the overall complex
        pdb_a_path = os.path.join(HOMs_directory, f"{pdb_id_a}.pdb")
        pdb_b_path = os.path.join(HOMs_directory, f"{pdb_id_b}.pdb")
        chain_a_length = len(load_sequence(pdb_a_path))
        chain_b_length = len(load_sequence(pdb_b_path))
        # Percent identity
        nw_score, percent_identity = calculate_percent_identity(f"{pdb_id_a}", pdb_a_path, f"{pdb_id_b}", pdb_b_path, os.path.join(temp_path, "Global_Alignment.txt"))
        # TM-Align Score
        structure_similarity = calculate_TMScore(pdb_a_path, pdb_b_path, os.path.join(temp_path, "TMAlign_Output.txt"), alignment=os.path.join(temp_path, "Global_Alignment.txt"))
        # Save the data
        line = f"{pdb_id_a},{pdb_id_b},{chain_a_length},{chain_b_length},{abs(chain_a_length-chain_b_length)},{nw_score},{nw_score/max(chain_a_length, chain_b_length)},{percent_identity},{structure_similarity}"
        print("    "+line)
        random_comparisons = random_comparisons.append({
            "PDB_ID_1":pdb_id_a,
            "PDB_ID_2":pdb_id_b,
            "PDB_1_Length":chain_a_length,
            "PDB_2_Length":chain_b_length,
            "Length_Difference":abs(chain_a_length - chain_b_length),
            "Alignment_Score":nw_score,
            "Alignment_Score_Adjusted":nw_score/max(chain_a_length, chain_b_length),
            "PID":percent_identity,
            "TM":structure_similarity}, ignore_index=True)


    random_comparisons.to_csv(output_file_path)



        



                    






