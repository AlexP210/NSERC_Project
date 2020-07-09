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
    for species_folder in os.listdir(root_directory):

        chains_directory = os.path.join(root_directory, species_folder, "NonRedundant_Chains")
        assemblies_directory = os.path.join(root_directory, species_folder, "Selected_Biological_Assemblies")
        temp_path = os.path.join(root_directory, species_folder, "_Temp")
        output_file_path = os.path.join(root_directory, species_folder, "RandomComparisons.csv")

        chain_files = os.listdir(chains_directory)
        chain_names = [cn.split(".")[0] for cn in chain_files]
        counter = 0
        pairings_examined = 0
        with open(output_file_path, "w") as output_file:
            output_file.write("PDB_ID_1,PDB_ID_2,Chain_1,Chain_2,Chain_1_Length,Chain_2_Length,Alignment_Score,PID,TM\n")
            for chain_idx_a in range(len(chain_names)):
                for chain_idx_b in range(chain_idx_a):

                    if pairings_examined >= n_pairings: 
                        sys.exit()

                    pdb_id_a = chain_names[chain_idx_a][:4].lower()
                    chain_letter_a = chain_names[chain_idx_a][5].upper()
                    chain_a_filename = f"{pdb_id_a}_{chain_letter_a}.ent"
                    pdb_a_filename = f"{pdb_id_a}.cif"

                    pdb_id_b = chain_names[chain_idx_b][:4].lower()
                    chain_letter_b = chain_names[chain_idx_b][5].upper()
                    chain_b_filename = f"{pdb_id_b}_{chain_letter_b}.ent"
                    pdb_b_filename = f"{pdb_id_b}.cif"

                    if counter%1 == 0: # If the files exist, and this is the hundredth pairing examined, then analyze
                        pairings_examined += 1
                        print(f"Processing: {chain_a_filename} vs {chain_b_filename} # {pairings_examined}")
                        # Get the paths to the chains, and the overall complex
                        chain_a_path = os.path.join(chains_directory, f"{pdb_id_a}_{chain_letter_a}.ent")
                        chain_b_path = os.path.join(chains_directory, f"{pdb_id_b}_{chain_letter_b}.ent")
                        # Percent identity
                        nw_score, percent_identity = calculate_percent_identity(f"{pdb_id_a}_{chain_letter_a}", chain_a_path, f"{pdb_id_b}_{chain_letter_b}", chain_b_path, os.path.join(temp_path, "Global_Alignment.txt"))
                        # TM-Align Score
                        structure_similarity = calculate_TMScore(chain_a_path, chain_b_path, os.path.join(temp_path, "TMAlign_Output.txt"), alignment=os.path.join(temp_path, "Global_Alignment.txt"))
                        # Save the data
                        line = f"{pdb_id_a},{pdb_id_b},{chain_letter_a},{chain_letter_b},{len(load_sequence(chain_a_path))},{len(load_sequence(chain_b_path))},{nw_score},{percent_identity},{structure_similarity}"
                        print("    "+line)
                        output_file.write(line + "\n")

                    
                    # Increment our counter
                    counter += 1


                    






