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
import My_Library as ml

from My_Library import *

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) < 3:
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <Species Root Directory> <Symmetry Groups>\n")
        sys.exit()
    root_directory = sys.argv[1]
    symmetry_groups = sys.argv[2:]
    print("1. COMPARING CHAINS")
    for species_fol in os.listdir(root_directory):
        print(species_fol)
        # Set all relevant folders for processing this species' data
        temp_dir = os.path.join(root_directory, species_fol, "_Temp")
        if not os.path.exists(temp_dir): os.mkdir(temp_dir)
        HOMs_directory = os.path.join(root_directory, species_fol, "Clean_Structures")
        assemblies_directory = os.path.join(root_directory, species_fol, "Selected_Biological_Assemblies")
        csv_path = os.path.join(root_directory, species_fol, "HOMs.csv")

        with open(csv_path, "w") as output_file:
            # Create and write the header for the csv
            header = "PDB_ID_1,PDB_ID_2,PDB_1_Length,PDB_2_Length,Length_Difference,Alignment_Score,Alignment_Score_Adjusted,PID,TM,TM_Rotation,TM_Translation"
            # for symmetry_group in symmetry_groups:
            #     header += f", {symmetry_group.upper()}_RMSD"
            # header += "\n"
            output_file.write(header)


            # id_to_chains = {}
            # for filename in os.listdir(chains_directory):
            #     pdb_id = filename.split("_")[0]
            #     chain = filename.split("_")[-1].split(".")[0]
            #     if pdb_id not in id_to_chains:
            #         id_to_chains[pdb_id] = [chain]
            #     else:
            #         id_to_chains[pdb_id].append(chain)

            pdb_file_list = [os.path.join(HOMs_directory, i) for i in os.listdir(HOMs_directory)]
            for pdb_file_1_idx in range(len(pdb_file_list)):
                for pdb_file_2_idx in range(pdb_file_1_idx, len(pdb_file_list)):
                    pdb_1_name = os.path.basename(pdb_file_list[pdb_file_1_idx]).split(".")[0]
                    pdb_1_path = pdb_file_list[pdb_file_1_idx]
                    pdb_2_name = os.path.basename(pdb_file_list[pdb_file_2_idx]).split(".")[0]
                    pdb_2_path = pdb_file_list[pdb_file_2_idx]
                    nw_score, percent_identity = calculate_percent_identity(pdb_1_name, pdb_1_path, pdb_2_name, pdb_2_path, os.path.join(temp_dir, "Global_Alignment.fasta"))
                    structure_similarity = calculate_TMScore(pdb_1_path, pdb_2_path, os.path.join(temp_dir, "TMAlign_Output.txt"), alignment = os.path.join(temp_dir, "Global_Alignment.fasta"), matrix_out=os.path.join(temp_dir, "TMAlign_Matrix.txt"))
                    x_rot, y_rot, z_rot, rotation_angle, x_trans, y_trans, z_trans, trans_mag = ml.calculate_rotation_angles(os.path.join(temp_dir, "TMAlign_Matrix.txt"))
                    length_a = len(load_sequence(pdb_1_path))
                    length_b = len(load_sequence(pdb_2_path))
                    line = f"{pdb_1_name},{pdb_2_name},{length_a},{length_b},{abs(length_a - length_b)},{nw_score},{nw_score/max(length_a, length_b)},{percent_identity},{structure_similarity},{rotation_angle},{trans_mag}"
                    print(line)
                    output_file.write(line + "\n")
