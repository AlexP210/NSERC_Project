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
    if len(sys.argv) < 3:
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <Species Root Directory> <Symmetry Groups>\n")
        sys.exit()
    root_directory = sys.argv[1]
    symmetry_groups = sys.argv[2:]
    print("1. COMPARING CHAINS")
    for species_fol in os.listdir(root_directory)[:1]:
        print(species_fol)
        # Set all relevant folders for processing this species' data
        temp_dir = os.path.join(root_directory, species_fol, "_Temp")
        if not os.path.exists(temp_dir): os.mkdir(temp_dir)
        interface_directory = os.path.join(root_directory, species_fol, "Separated_Interfaces")
        assemblies_directory = os.path.join(root_directory, species_fol, "Selected_Biological_Assemblies")
        csv_path = os.path.join(root_directory, species_fol, "Interfaces.csv")

        with open(csv_path, "w") as output_file:
            # Create and write the header for the csv
            header = "PDB_ID,Interface_1,Interface_2,Interface_1_Length,Interface_2_Length,Alignment_Score,PID,TM"
            for symmetry_group in symmetry_groups:
                header += f", {symmetry_group.upper()}_RMSD"
            header += "\n"
            output_file.write(header)


            id_to_chains = {}
            for filename in os.listdir(interface_directory):
                interface_id = filename.split("_")[0][:6]
                chain = filename.split("_")[-1].split(".")[0]
                if interface_id not in id_to_chains:
                    id_to_chains[interface_id] = [chain]
                else:
                    id_to_chains[interface_id].append(chain)

            for pdb_id, chain_letters in id_to_chains.items():
                for sequence_a_idx in range(len(chain_letters)-1):
                    for sequence_b_idx in range(sequence_a_idx+1, len(chain_letters)):
                        chain_a_letter = chain_letters[sequence_a_idx]
                        chain_b_letter = chain_letters[sequence_b_idx]
                        # Get the paths to the chains, and the overall complex
                        chain_a_path = os.path.join(interface_directory, f"{pdb_id}_{chain_a_letter}.ent")
                        chain_b_path = os.path.join(interface_directory, f"{pdb_id}_{chain_b_letter}.ent")
                        # complex_path = os.path.join(assemblies_directory, f"{pdb_id}.cif") - WORK
                        # Percent identity
                        nw_score, percent_identity = calculate_percent_identity(f"{pdb_id}_{chain_a_letter}", chain_a_path, f"{pdb_id}_{chain_b_letter}", chain_b_path, os.path.join(temp_dir, "Global_Alignment.fasta"))
                        # TM-Align Score
                        structure_similarity = calculate_TMScore(chain_a_path, chain_b_path, os.path.join(temp_dir, "TMAlign_Output.txt"), alignment = os.path.join(temp_dir, "Global_Alignment.fasta"))
                        # Symmetry RMSD
                        # symmetry_rmsds = calculate_symmetry_rmsd(complex_path, symmetry_groups, os.path.join(temp_dir, "AnAnaS_Output.txt")) - WORK
                        # Save the data
                        line = f"{pdb_id[0:4]},{chain_letters[sequence_a_idx]},{chain_letters[sequence_b_idx]},{len(load_sequence(chain_a_path))},{len(load_sequence(chain_b_path))},{nw_score},{percent_identity},{structure_similarity}"
                        # for symmetry_rmsd in symmetry_rmsds:
                        #     line += f", {symmetry_rmsd}"
                        print("    "+line)
                        output_file.write(line + "\n")

                    


    




