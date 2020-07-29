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

    # # for record in SeqIO.parse( os.path.join(directory, f"{pdb_id}.cif"), "cif-seqres" ):
    # #     print(record.name)
    # #     if record.id.split(":")[-1] == chain:
    # #         return record.seq

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(f"{pdb_id}_{chain}", os.path.join(directory, f"{pdb_id}.cif"))
    sequence = ""
    for chain_structure in structure.get_chains():
        if chain_structure.get_id() == chain:
            for residue in chain_structure.get_residues():
                sequence += one_letter_code(residue.get_resname())
            return sequence
        

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
        translation_vector = []
        for row_text in matrix_rows_text:
            matrix_row = []
            for entry in row_text:
                if entry != "": matrix_row.append(float(entry))
            matrix.append(matrix_row[1:])
            translation_vector.append(matrix_row[0])
    return matrix, translation_vector

def calculate_rotation_angles(path_to_matrix):
    # Parse the matrix
    rotation_matrix, translation_vector = _parse_rotation_matrix(path_to_matrix)
    # Make sure it's right
    det = np.linalg.det(rotation_matrix)
    assert abs(det - 1) < 0.001, f"Found rotation matrix with determinant {det}"
    rotation = Rotation.from_matrix(rotation_matrix)
    # Unpack the angles
    x_rot, y_rot, z_rot = rotation.as_euler("xyz") # Euler angle representation
    rotation_angle = rotation.magnitude() # Overall rotation angle
    for angle in (x_rot, y_rot, z_rot, rotation_angle):
        angle = angle*(180/math.pi)
    # Unpack the translation vector
    x_trans, y_trans, z_trans = translation_vector[0], translation_vector[1], translation_vector[2]
    trans_mag = ( x_trans**2 + y_trans**2 + z_trans**2 )**(0.5)
    return x_rot, y_rot, z_rot, rotation_angle, x_trans, y_trans, z_trans, trans_mag

def calculate_c2Rmsd(path_to_ananas):
    with open(path_to_ananas, "r") as ananas_output:
        for line in ananas_output:
            if line[:13] == "No symmetries":
                return -1
            if line[:12] == "Average RMSD":
                return float(line.split(" : ")[-1])


if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) != 6:
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <Chains Directory> <Biological Assemblies Directory> <CIFs Directory> <output.csv> <Temp Directory>\n")
        sys.exit()
    chains_directory = sys.argv[1]
    assemblies_directory = sys.argv[2]
    cif_directory = sys.argv[3]
    csv_path = sys.argv[4]
    temp_dir = sys.argv[5]
    if not os.path.exists(temp_dir): os.mkdir(temp_dir)

    # Open the output file and write the sequence similarities
    with open(csv_path, "w") as output_file:
        output_file.write("PDB_ID, Chain_1, Chain_2, Chain_1_Length, Chain_2_Length, Percent_Identity_of_Global_Alignment, TM_Score, Identity_Of_Structural_Alignment, X_Rotation, Y_Rotation, Z_Rotation, Total_Rotation, X_Translation, Y_Translation, Z_Translation, Total_Translation, C2_RMSD\n")
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
            chain_sequences = [ load_sequence(pdb_id, chain_letter, assemblies_directory) for chain_letter in chain_letters ]
            for sequence_a_idx in range(len(chain_sequences)-1):
                for sequence_b_idx in range(sequence_a_idx+1, len(chain_sequences)):
                    # Calculate Sequence alignment
                    sequence_a = chain_sequences[sequence_a_idx]
                    sequence_b = chain_sequences[sequence_b_idx]
                    chain_a_letter = chain_letters[sequence_a_idx]
                    chain_b_letter = chain_letters[sequence_b_idx]
                    seq_score = global_alignment_score(sequence_a, sequence_b)
                    chain_a_path = os.path.join(chains_directory, f"{pdb_id}_{chain_a_letter}.ent")
                    chain_b_path = os.path.join(chains_directory, f"{pdb_id}_{chain_b_letter}.ent")
                    output_path = csv_path
                    matrix_path = os.path.join(temp_dir, "matrix.txt")
                    ananas_output_path = os.path.join(temp_dir, "ananas.txt")
                    tm_align_path = "TMalign"
                    ananas_path = "ananas"
                    output_path = os.path.join(temp_dir, "struc_alignment.txt")
                    cif_filename = os.path.join(assemblies_directory, f"{pdb_id}.cif")
                    call = f'{tm_align_path} {chain_a_path} {chain_b_path} -m {matrix_path} > {output_path} & {ananas_path} {cif_filename} -C 100 c3 > {ananas_output_path}'
                    os.system(call)
                    # Calculate Structure Alignment - extract the TM score from the output, and remove it
                    with open(output_path, "r") as alignment_output:
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
                    x_rot, y_rot, z_rot, rotation_angle, x_trans, y_trans, z_trans, trans_mag = calculate_rotation_angles(os.path.join(temp_dir, "matrix.txt"))
                    # Calculate the RMSD for C2 alignment
                    c2_alignment = calculate_c2Rmsd(os.path.join(temp_dir, "ananas.txt"))
                    print(f"{pdb_id}, {chain_letters[sequence_a_idx]}, {chain_letters[sequence_b_idx]}, {len(sequence_a)}, {len(sequence_b)}, {seq_score}, {TM_score}, {identity}, {x_rot}, {y_rot}, {z_rot}, {rotation_angle}, {x_trans}, {y_trans}, {z_trans}, {trans_mag}, {c2_alignment}")
                    output_file.write(f"{pdb_id}, {chain_letters[sequence_a_idx]}, {chain_letters[sequence_b_idx]}, {len(sequence_a)}, {len(sequence_b)}, {seq_score}, {TM_score}, {identity}, {x_rot}, {y_rot}, {z_rot}, {rotation_angle}, {x_trans}, {y_trans}, {z_trans}, {trans_mag}, {c2_alignment}\n")

                    


    




