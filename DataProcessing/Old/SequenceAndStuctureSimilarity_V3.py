# Changes over V2:
# 1. Also calcualates the rotation matrix used to align the structures, and records rotation angles

from Bio.PDB.PDBParser import PDBParser
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

def one_letter_code(residue_name):
    conversion = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    
    if residue_name in conversion:
        return conversion[residue_name]
    else:
        return ""

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

def _parse_rotation_matrix(path_to_matrix):
    with open(path_to_matrix, "r") as matrix_text:
        text_lines = matrix_text.readlines()
        matrix_rows_text = [text_lines[2:5][i][1:-1].split(" ")[1:] for i in range(3)]
        matrix = []
        for row_text in matrix_rows_text:
            matrix_row = []
            for entry in row_text:
                if entry != "": matrix_row.append(float(entry))
            matrix.append(matrix_row[1:])
    return matrix

def calculate_rotation_angles(path_to_matrix):
    rotation_matrix = _parse_rotation_matrix("C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\Stuctural_Similarity\\matrix.txt")
    det = np.linalg.det(rotation_matrix)
    assert abs(det - 1) < 0.001, f"Found rotation matrix with determinant {det}"
    rotation = Rotation.from_matrix(rotation_matrix)
    z, y, x = rotation.as_euler("zyx") # Euler angle representation
    rotation_angle = rotation.magnitude() # Overall rotation angle
    for angle in (z, y, x, rotation_angle):
        angle = angle*180/(math.pi)
    return z, y, x, rotation_angle

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) != 5:
        print(f"Usage: $ python SequenceAndStructureSimilarity_V3.py <PDB dir> <chains dir> <chains.txt> <output.txt>")
        sys.exit()
    pdb_directory = sys.argv[1]
    chains_directory = sys.argv[2]
    chain_names_path = sys.argv[3]
    csv_path = sys.argv[4]

    # First read in the chain names from the PDBs.txt file, 
    # create a dictionary with (key, val) = (PDB ID, Chains of Interest)
    with open(chain_names_path, "r") as chain_names:
        chain_dict = {}
        for unparsed_chain_name in chain_names:
            pdb_id = unparsed_chain_name[:4].lower()
            chain_letter = unparsed_chain_name[4].upper()
            if pdb_id not in chain_dict:
                chain_dict[pdb_id] = [chain_letter]
            else:
                chain_dict[pdb_id].append(chain_letter)

    # Open the output file and write the sequence similarities
    with open(csv_path, "w") as output_file:
        output_file.write("PDB_ID, Chain_1, Chain_2, Percent_Identity_of_Global_Alignment, TM_Score, Identity_Of_Structural_Alignment, Z_Rotation, Y_Rotation, X_Rotation, Total_Rotation\n")
        for pdb_id, chain_letters in chain_dict.items():
            if not os.path.isfile(os.path.join(pdb_directory, f"pdb{pdb_id}.ent")):
                print(f"WARNING: pdb{pdb_id}.ent not found. Check downloads.")
            else: # The PDB was downloaded properly
                chain_sequences = [ load_sequence(pdb_id, chain_letter, pdb_directory) for chain_letter in chain_letters ]
                for sequence_a_idx in range(len(chain_sequences)-1):
                    for sequence_b_idx in range(sequence_a_idx+1, len(chain_sequences)):
                        # Calculate Sequence alignment
                        sequence_a = chain_sequences[sequence_a_idx]
                        sequence_b = chain_sequences[sequence_b_idx]
                        seq_score = global_alignment_score(sequence_a, sequence_b)
                        # Calculate structure alignment - Execute TMalign program through WSL, send outputs to current directory
                        # # Paths to relevant files in WSL
                        chain_a_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/pdb{pdb_id}_{chain_letters[sequence_a_idx]}.ent"
                        chain_b_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/pdb{pdb_id}_{chain_letters[sequence_b_idx]}.ent"
                        output_path = f"../../mnt/c/Users/alexp/Desktop/School/3rd_Year/NSERC/StructuralSequenceSimilarity/Scripts/Stuctural_Similarity/struc_alignment.txt"
                        matrix_path = f"../../mnt/c/Users/alexp/Desktop/School/3rd_Year/NSERC/StructuralSequenceSimilarity/Scripts/Stuctural_Similarity/matrix.txt"
                        tm_align_path = f"./Protein_Complexes/TMAlign/TMalign"
                        call = f'ubuntu run "cd ~;{tm_align_path} {chain_a_path} {chain_b_path} -m {matrix_path} > {output_path}"'
                        os.system(call)
                        # Calculate Structure Alignment - extract the TM score from the output, and remove it
                        with open('C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\Stuctural_Similarity\\struc_alignment.txt', "r") as alignment_output:
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
                        z, y, x, rotation_angle = calculate_rotation_angles("C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\Stuctural_Similarity\\matrix.txt")
                        print(f"{pdb_id}, {chain_letters[sequence_a_idx]}, {chain_letters[sequence_b_idx]}, {seq_score}, {TM_score}, {identity}, {z}, {y}, {x}, {rotation_angle}")
                        output_file.write(f"{pdb_id}, {chain_letters[sequence_a_idx]}, {chain_letters[sequence_b_idx]}, {seq_score}, {TM_score}, {identity}, {z}, {y}, {x}, {rotation_angle}\n")
    # Remove the temporary files used for TMAlign output
    os.remove('C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\Stuctural_Similarity\\struc_alignment.txt')
    os.remove("C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\Stuctural_Similarity\\matrix.txt")

                    


    




