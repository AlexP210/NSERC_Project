import Sequence_Identities_NW_V2 as si
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
import os.path
import os
import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import align, format_alignment
from SequenceAndStuctureSimilarity_V3 import calculate_rotation_angles

if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) < 6:
        print(f"Usage: $ python {__file__} <chains.txt> <pdb directory> <chains directory> <output.csv> <Number Pairings>")
        sys.exit()
    chains_list_path = sys.argv[1]
    pdb_directory = sys.argv[2]
    chains_directory = sys.argv[3]
    output_file_path = sys.argv[4]
    n_pairings = sys.argv[5]

    # Go over the list of chains, and compare any two chains that:
    #   1. Have both chains existing
    #   2. Were not previously compared
    #   3. Are not identical
    #   4. Have index as multiple of 500
    with open(chains_list_path, "r") as chains_list:
        with open(output_file_path, "w") as output_file:
            output_file.write("PDB_ID_1, PDB_ID_2, Chain_1, Chain_2, Percent_Identity_of_Global_Alignment, TM_Score, Identity_Of_Structural_Alignment, z, y, x, Total_Rotation\n")
            chain_names = chains_list.readlines()
            counter = 0
            pairings_examined = 0
            for chain_idx_a in range(len(chain_names)):
                for chain_idx_b in range(chain_idx_a):

                    if pairings_examined >= n_pairings: sys.exit()

                    pdb_id_a = chain_names[chain_idx_a][:4].lower()
                    chain_letter_a = chain_names[chain_idx_a][4].upper()
                    chain_a_filename = f"pdb{pdb_id_a}_{chain_letter_a}.ent"
                    pdb_a_filename = f"pdb{pdb_id_a}.ent"

                    pdb_id_b = chain_names[chain_idx_b][:4].lower()
                    chain_letter_b = chain_names[chain_idx_b][4].upper()
                    chain_b_filename = f"pdb{pdb_id_b}_{chain_letter_b}.ent"
                    pdb_b_filename = f"pdb{pdb_id_b}.ent"

                    if not (os.path.isfile( os.path.join(pdb_directory, pdb_a_filename) ) and os.path.isfile( os.path.join(pdb_directory, pdb_b_filename) )):
                        print(f"Warning: {chain_a_filename} or {chain_b_filename} nonexistent.")
                    elif counter%500 == 0: # If the files exist, and this is the hundredth pairing examined, then analyze
                        # Get percent identity
                        chain_a_sequence = si.load_sequence(pdb_id_a, chain_letter_a, pdb_directory)
                        chain_b_sequence = si.load_sequence(pdb_id_b, chain_letter_b, pdb_directory)
                        seq_score = si.global_alignment_score(chain_a_sequence, chain_b_sequence)
                        # Get Structural Alignment Score
                        chain_a_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/{chain_a_filename}" # Need to fix these paths - If call script from WSL, don't need to worry about Windows/Linux path conversions. Could also bring TMalign into windows
                        chain_b_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/{chain_b_filename}"
                        temp_path = f"../../mnt/c/Users/alexp/Desktop/School/3rd_Year/NSERC/StructuralSequenceSimilarity/Scripts/Stuctural_Similarity/struc_alignment.txt"
                        tm_align_path = f"./Protein_Complexes/TMAlign/TMalign"
                        call = f'ubuntu run "cd ~;{tm_align_path} {chain_a_path} {chain_b_path}>{temp_path}"'
                        os.system(call)
                        with open('C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\Stuctural_Similarity\\struc_alignment.txt', "r") as alignment_output:
                            TM_scores = []
                            for line in alignment_output: # Parse the TM-Align output
                                if len(line) > 7: # Avoid index errors with short lines in output
                                    if line[0:7] == "Aligned":
                                        all_info = line.split(", ")
                                        identity_data = all_info[2].split("=")
                                        identity = float(identity_data[-1][1:])
                                    if line[0:8] == "TM-score":
                                        TM_score = float(line.split(" ")[1])
                                        TM_scores.append(TM_score)
                        z, y, x, total_rotation = calculate_rotation_angles('C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\Stuctural_Similarity\\struc_alignment.txt')
                        os.remove('C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\Stuctural_Similarity\\struc_alignment.txt')
                        TM_score = min(TM_scores)
                        pairings_examined += 1
                        # Write the output
                        print(f"Pairing #{pairings_examined}: {pdb_id_a}, {pdb_id_b}, {chain_letter_a}, {chain_letter_b}, {seq_score}, {TM_score}, {identity}, {z}, {x}, {y}, {total_rotation}")
                        output_file.write(f"{pdb_id_a}, {pdb_id_b}, {chain_letter_a}, {chain_letter_b}, {seq_score}, {TM_score}, {identity}, {z}, {x}, {y}, {total_rotation}\n")
                    
                    # Increment our counter
                    counter += 1


                    






