# Changes over V4:
# 1. Does not need to use WSL

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
import os.path
import os
import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import align, format_alignment
import numpy as np
from math import atan2, acos
from scipy.spatial.transform import Rotation
import math
import time

from My_Library import *

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) != 5:
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <Chains Directory> <Biological Assemblies Directory> <output.csv> <Temp Directory>\n")
        sys.exit()
    chains_directory = sys.argv[1]
    assemblies_directory = sys.argv[2]
    csv_path = sys.argv[3]
    temp_dir = sys.argv[4]
    if not os.path.exists(temp_dir): os.mkdir(temp_dir)

    # Open the output file and write the sequence similarities
    with open(csv_path, "w") as output_file:
        output_file.write("PDB_ID, Chain_1, Chain_2, Chain_1_Length, Chain_2_Length, Percent_Identity_of_Global_Alignment, TM_Score, Identity_Of_Structural_Alignment, X_Rotation, Y_Rotation, Z_Rotation, Total_Rotation, X_Translation, Y_Translation, Z_Translation, Total_Translation, C2_RMSD, D2_RMSD\n")
        id_to_chains = {}
        for filename in os.listdir(chains_directory):
            pdb_id = filename.split("_")[0]
            chain = filename.split("_")[-1].split(".")[0]
            if pdb_id not in id_to_chains:
                id_to_chains[pdb_id] = [chain]
            else:
                id_to_chains[pdb_id].append(chain)

        for pdb_id, chain_letters in id_to_chains.items():
            print(f"\nProcessing {pdb_id}...")
            for sequence_a_idx in range(len(chain_letters)-1):
                for sequence_b_idx in range(sequence_a_idx+1, len(chain_letters)):
                    chain_a_letter = chain_letters[sequence_a_idx]
                    chain_b_letter = chain_letters[sequence_b_idx]
                    # Get the paths to the chains, and the overall complex
                    chain_a_path = os.path.join(chains_directory, f"{pdb_id}_{chain_a_letter}.ent")
                    chain_b_path = os.path.join(chains_directory, f"{pdb_id}_{chain_b_letter}.ent")
                    complex_path = os.path.join(assemblies_directory, f"{pdb_id}.cif")
                    # Percent identity
                    percent_identity = calculate_percent_identity(chain_a_path, chain_b_path)
                    # TM-Align Score
                    structure_similarity = calculate_TMScore(chain_a_path, chain_b_path)
                    # Symmetry RMSD
                    c2_rmsd, d2_rmsd = calculate_symmetry_rmsd(complex_path, ("c2", "d2"))
                    # Save the data
                    print(f"{pdb_id}, {chain_letters[sequence_a_idx]}, {chain_letters[sequence_b_idx]}, {len(load_sequence(chain_a_path))}, {len(load_sequence(chain_b_path))}, {percent_identity}, {percent_identity}, {c2_rmsd}, {d2_rmsd}")
                    output_file.write(f"{pdb_id}, {chain_letters[sequence_a_idx]}, {chain_letters[sequence_b_idx]}, {len(load_sequence(chain_a_path))}, {len(load_sequence(chain_b_path))}, {percent_identity}, {percent_identity}, {c2_rmsd}, {d2_rmsd}\n")

                    


    




