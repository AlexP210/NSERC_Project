# Changes over V1:
#   1. calculates the rotation matrix
#   2. Excludes pairings corresponding to real heteromers
#   3. Does Ananas

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
import os.path
import os
import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import align, format_alignment
from HeteromerComparisons_V4 import calculate_rotation_angles
from HeteromerComparisons_V4 import calculate_c2Rmsd

def load_sequence(pdb_id, chain, directory):
    for record in SeqIO.parse( os.path.join(directory, f"pdb{pdb_id.lower()}.ent"), "pdb-seqres" ):
        if record.id.split(":")[-1].upper() == chain.upper():
            return record.seq

def global_alignment_score(sequence_a, sequence_b):
    alignment = align.globaldx(sequence_a, sequence_b, matlist.blosum62, one_alignment_only=True)[0]
    seq_score = 0
    for i in range(len(alignment[0])):
        if alignment[0][i] == alignment[1][i] and alignment[0][i] != "-":
            seq_score += 1
    seq_score /= max(len(sequence_a), len(sequence_b))
    return seq_score

if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) < 6:
        print(f"Usage: $ python {__file__.split('/')[-1]} <chains.txt> <pdb directory> <chains directory> <output.csv> <Number Pairings>")
        sys.exit()
    chains_list_path = sys.argv[1]
    pdb_directory = sys.argv[2]
    chains_directory = sys.argv[3]
    output_file_path = sys.argv[4]
    n_pairings = int(sys.argv[5])

    # Go over the list of chains, and compare any two chains that:
    #   1. Have both chains existing
    #   2. Were not previously compared
    #   3. Are not identical
    #   4. Have index as multiple of 500
    with open(chains_list_path, "r") as chains_list:
        with open(output_file_path, "w") as output_file:
            output_file.write("PDB_ID_1, PDB_ID_2, Chain_1, Chain_2, Chain_1_Length, Chain_2_Length, Percent_Identity_of_Global_Alignment, TM_Score, Identity_Of_Structural_Alignment, X_Rotation, Y_Rotation, Z_Rotation, Total_Rotation, X_Translation, Y_Translation, Z_Translation, Total_Translation, C2_RMSD\n")
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
                        pass
                    elif counter%500 == 0: # If the files exist, and this is the hundredth pairing examined, then analyze
                        pairings_examined += 1
                        # Get percent identity
                        chain_a_sequence = load_sequence(pdb_id_a, chain_letter_a, pdb_directory)
                        chain_b_sequence = load_sequence(pdb_id_b, chain_letter_b, pdb_directory)
                        seq_score = global_alignment_score(chain_a_sequence, chain_b_sequence)
                        # Get Structural Alignment Score
                        chain_a_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/{chain_a_filename}" # Need to fix these paths - If call script from WSL, don't need to worry about Windows/Linux path conversions. Could also bring TMalign into windows
                        chain_b_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/{chain_b_filename}"
                        temp_path = f"../../mnt/c/Users/alexp/Desktop/School/3rd_Year/NSERC/StructuralSequenceSimilarity/Scripts/Stuctural_Similarity/struc_alignment.txt"
                        tm_align_path = f"./Protein_Complexes/TMAlign/TMalign"
                        ananas_path = f"./Protein_Complexes/AnAnaS/ananas"
                        call = f'ubuntu run "cd ~;{tm_align_path} {chain_a_path} {chain_b_path}>{temp_path}; "'
                        os.system(call)
                         # Calculate Structure Alignment - extract the TM score from the output, and remove it
                        with open('C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\DataProcessing\\Temp\\struc_alignment.txt', "r") as alignment_output:
                            TM_scores = []
                            for line in alignment_output:
                                if line[0:2] == "Al":
                                    all_info = line.split(", ")
                                    identity_data = all_info[2].split("=")
                                    identity = float(identity_data[-1][1:])
                                if line[0:2] == "TM":
                                    TM_score_data = line.split("= ")
                                    TM_score = float(TM_score_data[1].split(" ")[0])
                                    TM_scores.append(TM_score)
                        TM_score = min(TM_scores)
                        # Calculate structure alignment - extract the rotation angles from the matrix output, and remove it
                        x_rot, y_rot, z_rot, rotation_angle, x_trans, y_trans, z_trans, trans_mag = calculate_rotation_angles("C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\DataProcessing\\Temp\\matrix.txt")
                        # Calculate the RMSD for C2 alignment
                        c2_alignment = calculate_c2Rmsd("C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\DataProcessing\\Temp\\ananas.txt")
                        print(f"{pdb_id_a}, {pdb_id_b}, {chain_letter_a}, {chain_letter_b}, {len(chain_a_sequence)}, {len(chain_b_sequence)}, {seq_score}, {TM_score}, {identity}, {x_rot}, {y_rot}, {z_rot}, {rotation_angle}, {x_trans}, {y_trans}, {z_trans}, {trans_mag}, {c2_alignment}")
                        output_file.write(f"{pdb_id_a}, {pdb_id_b}, {chain_letter_a}, {chain_letter_b}, {len(chain_a_sequence)}, {len(chain_b_sequence)}, {seq_score}, {TM_score}, {identity}, {x_rot}, {y_rot}, {z_rot}, {rotation_angle}, {x_trans}, {y_trans}, {z_trans}, {trans_mag}, {c2_alignment}\n")
                    
                    # Increment our counter
                    counter += 1


                    






