from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
from Bio.PDB import MMCIFIO, PDBIO
from Bio.PDB import Select
import os.path
import os
import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Align import substitution_matrices
from Bio.pairwise2 import align, format_alignment
import numpy as np
from math import atan2, acos
from scipy.spatial.transform import Rotation
import math
import time
import math as m
import numpy as np
import scipy


def calculate_symmetry_rmsd(pdb_filepath, symmetry_groups, temp_path):
    """Calculates RMSD for a specific symmetry group, using AnAnaS.exe

    :param pdb_filepath: Path of the PDB structure to analyze
    :type pdb_filepath: String
    :param symmetry_groups: Set of symmetry groups to use
    :type symmetry_groups: Iterable(String)
    :param temp_path: Path to temporary .txt file to output, defaults to "Ananas_Output.txt"
    :type temp_path: str, optional
    :return: RMSDs for each symmetry group specified
    :rtype: Iterable(float)
    """
    call = f"ananas {pdb_filepath} -C 100 "
    for sym in symmetry_groups:
        call += sym
    call += f"> {temp_path}"
    os.system(call)
    next_rmsd = False
    rmsds = {}
    with open(temp_path, "r") as ananas_output:
        for line in ananas_output:
            if line[:13] == "No symmetries":
                return [-1]
            if line.split(" : ")[0] == "Symmetry group" and line.split(" : ")[-1][:2] in symmetry_groups:
                next_rmsd = True
                next_group = line.split(" : ")[-1][:2]
            if line[:12] == "Average RMSD" and next_rmsd:
                rmsds[next_group] = float(line.split(" : ")[-1])
    os.remove(temp_path)
    return [rmsds[group] for group in symmetry_groups]

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

def global_alignment_score(sequence_a, sequence_b):
    alignment = align.globaldx(sequence_a, sequence_b, matlist.pam60, one_alignment_only=True)[0]
    seq_score = 0
    for i in range(len(alignment[0])):
        if alignment[0][i] == alignment[1][i] and alignment[0][i] != "-":
            seq_score += 1
    seq_score /= max(len(sequence_a), len(sequence_b))
    return seq_score

def load_sequence(pdb_filepath):
    if pdb_filepath.split(".")[-1] == "cif": parser = MMCIFParser(QUIET=True)
    if pdb_filepath.split(".")[-1] in ("ent", "pdb"): parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_filepath.split("/")[-1], pdb_filepath)
    sequence = ""
    for chain_structure in structure.get_chains():
        for residue in chain_structure.get_residues():
            sequence += one_letter_code(residue.get_resname())
    return sequence

def one_letter_code(residue_name):
    conversion = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    
    if residue_name in conversion:
        return conversion[residue_name]
    else:
        return ""

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

def calculate_TMScore(pdb_1_path, pdb_2_path, temp_path, alignment = None):
    """Calculates the optimal structure alignment between two structures, and outputs the TM-Align score.

    :param pdb_1_path: Path to one structure
    :type pdb_1_path: String
    :param pdb_2_path: Path to other structure
    :type pdb_2_path: String
    :param temp_path: Path to temporary .txt to output, defaults to "TMalign_Output.txt"
    :type temp_path: str, optional
    :param alignment: Path to .fst file to use for the alignment
    :type alignment: str, optional
    :return: TM-Align score for optimal structural alignment
    :rtype: String
    """
    call = f"TMalign {pdb_1_path} {pdb_2_path}"
    if alignment != None:
        call += f" -I {alignment}"
    call += f" > {temp_path}"
    os.system(call)
    with open(temp_path, "r") as alignment_output:
        TM_scores = []
        for line in alignment_output:
            if line[0:2] == "Al" or line[0:5] == "ERROR":
                # This is the case that the structure alignment failed (chain alignment is too small)
                pass
            if line[0:2] == "TM":
                TM_score_data = line.split("= ")
                TM_score = float(TM_score_data[1].split(" ")[0])
                TM_scores.append(TM_score)
    if len(TM_scores) != 0:
        TM_score = min(TM_scores)
    else:
        TM_score = -1
    os.remove(temp_path)
    return TM_score

def get_chain_ids(path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(path.split("/")[-1].split("-")[0], path)
    chain_letters = []
    for chain in structure.get_chains():
        chain_letters.append(chain.get_id())
    return chain_letters

def calculate_percent_identity(pdb_1_name, pdb_1_path, pdb_2_name, pdb_2_path, save_alignment_path = None):
    """Calculates % identity for the global alignment between 2 protein sequences

    :param pdb_1_path: Path to structure file.
    :type pdb_1_path: String
    :param pdb_2_path: Path to other structure.
    :type pdb_2_path: String
    :return: (Alignment Score, Percent Identity)
    :rtype: Tuple<float>
    """
    # Create the aligner object
    pam70 = substitution_matrices.load("PAM70")
    aligner = Align.PairwiseAligner(substitution_matrix = pam70, open_gap_score=-5, extend_gap_score=-1)
    # Get the sequences
    sequence_a = load_sequence(pdb_1_path)
    sequence_b = load_sequence(pdb_2_path)
    # If the interface is too small, then return -1
    if len(sequence_a) == 0 or len(sequence_b) == 0:
        return -1, -1
    # Perform the alignment adn check the identity
    alignment = aligner.align(sequence_a, sequence_b)[0]
    aligned_subsequences = alignment.aligned
    nw_score = alignment.score
    percent_identity = 0
    for subsequence_idx in range(len(aligned_subsequences[0])):
        for residue_idx_a, residue_idx_b in zip(range(*aligned_subsequences[0][subsequence_idx]), range(*aligned_subsequences[1][subsequence_idx])):
            residue_a = sequence_a[residue_idx_a]
            residue_b = sequence_b[residue_idx_b]
            if residue_a == residue_b and not (residue_a == "-" or residue_b == "-"):
                percent_identity += 1
    percent_identity /= max(len(sequence_a), len(sequence_b))
    # Save the alignment as a fasta file if needed
    # To save the alignment to fasta, go through the alignment.aligned list, and write as many "-" as
    # there is space between the last residue in aligned subsequence n-1 and the first residue in aligned subsequence n.
    # Then, write the residues of aligned subsequence n
    if save_alignment_path != None:
        with open(save_alignment_path, "w") as fasta:
            sequence_a_aligned = ""
            sequence_b_aligned = ""
            last_subsequence_a = (0, 0)
            last_subsequence_b = (0, 0)
            for subsequence_idx in range(len(aligned_subsequences[0])):
                subsequence_a = aligned_subsequences[0][subsequence_idx]
                subsequence_b = aligned_subsequences[1][subsequence_idx]
                subsequence_a_skip = subsequence_a[0] - last_subsequence_a[1]
                subsequence_b_skip = subsequence_b[0] - last_subsequence_b[1]
                if subsequence_a_skip < subsequence_b_skip:
                    sequence_a_aligned += "-"*subsequence_b_skip + sequence_a[subsequence_a[0] : subsequence_a[1]]
                    sequence_b_aligned += sequence_b[last_subsequence_b[1] : subsequence_b[1]]
                elif subsequence_a_skip > subsequence_b_skip:
                    sequence_a_aligned += sequence_a[last_subsequence_a[1] : subsequence_a[1]]
                    sequence_b_aligned += "-"*subsequence_a_skip + sequence_b[subsequence_b[0] : subsequence_b[1]]
                else:
                    sequence_a_aligned += sequence_a[subsequence_a[0] : subsequence_a[1]]
                    sequence_b_aligned += sequence_b[subsequence_b[0] : subsequence_b[1]]
                last_subsequence_a = subsequence_a
                last_subsequence_b = subsequence_b
            # Write first sequence
            fasta.write(">"+pdb_1_name+"\n") 
            fasta.write(sequence_a_aligned)
            fasta.write("\n")
            # Write second sequence
            fasta.write(">"+pdb_2_name+"\n")
            fasta.write(sequence_b_aligned)
    # Return the identity, and score
    return nw_score, percent_identity

def gamma_probability(x, theta, k, mu):
    return 1/( m.factorial(k-1) * theta ** k ) * ((x - mu)**(k-1)) * m.e**(-(x-mu)/theta)

def log_gamma_probability(x, theta, k, mu):
    return -m.log(m.factorial(k-1)) - k*m.log(theta) + (k-1)*m.log(x-mu) - (x-mu)/theta

def fit_gamma(X):
    theta_hypotheses = [t/50 for t in range(1, 101)]
    k_hypotheses = [t/5 for t in range(1, 101)]
    mu_hypotheses = [i/100 for i in range(0, 101)]

    log_likelihood_space = np.ones((len(theta_hypotheses), len(k_hypotheses), len(mu_hypotheses)))
    max_log_likelihood = -float("inf")
    max_log_likelihood_idx = (None, None, None)
    for theta_idx in range(len(theta_hypotheses)):
        for k_idx in range(len(k_hypotheses)):
            for mu_idx in range(len(mu_hypotheses)):

                # Set the hypothesis
                theta = theta_hypotheses[theta_idx]
                k = k_hypotheses[k_idx]
                mu = mu_hypotheses[mu_idx]

                # Get the log-likelihood
                log_likelihood = 0
                for x in X:
                    try:
                        log_likelihood += m.log(gamma_probability(x, theta, k, mu))
                    except ValueError:
                        log_likelihood = -float("inf")
                if log_likelihood > max_log_likelihood: 
                    max_log_likelihood = log_likelihood
                    max_log_likelihood_idx = (theta_idx, k_idx, mu_idx)
                log_likelihood_space[theta_idx, k_idx, mu_idx] = log_likelihood
    best_theta, best_k, best_mu = max_log_likelihood_idx
    return lambda x: gamma_probability(x, theta_hypotheses[best_theta], best_k, mu_hypotheses[best_mu])

def p_value(X, prob):
    hist = {}
    for x in X:
        if not x in hist:
            hist[x] = 1
        else:
            hist[x] += 1
    sorted_hist = sorted(hist.items(), key=lambda i: i[0])
    total = 0
    for i in sorted_hist:
        total += i[1]
        if total / len(X) >= 0.95:
            return i[0]

def create_cdf(X):
    if len(X) > 0:
        return lambda x: sum([i <= x for i in X])/len(X)
    else: 
        return lambda x: 0

def get_interface(pdb_path_1, pdb_path_2, distance_cutoff, interface_path_1, interface_path_2):

    # Selector class
    class SelectInterface(Select):
        def __init__(self, other_structure, cutoff):
            self.other_structure = other_structure
            self.cutoff = cutoff
        def accept_residue(self, residue):
            for other_atom in self.other_structure.get_atoms():
                for atom in residue.get_atoms():
                    if np.linalg.norm(atom.get_coord() - other_atom.get_coord()) <= self.cutoff:
                        return 1
            return 0

    # Load in the structures, and create the Parser and IO objects to handle the input files
    parsers = []
    ios = []
    for input_file in (pdb_path_1, pdb_path_2):
        extension = os.path.basename(input_file).split(".")[-1]
        if extension == "cif": 
            parsers.append(MMCIFParser(QUIET=True))
            ios.append(MMCIFIO())
        elif extension == "ent": 
            parsers.append(PDBParser(QUIET=True))
            ios.append(PDBIO())
    structure_1 = parsers[0].get_structure("pdb_1", pdb_path_1)
    structure_2 = parsers[1].get_structure("pdb_2", pdb_path_2)

    # Create the selectors
    structure_1_selector = SelectInterface(structure_2, distance_cutoff)
    structure_2_selector = SelectInterface(structure_1, distance_cutoff)

    # Save the chains
    ios[0].set_structure(structure_1)
    ios[0].save(interface_path_1, select=structure_1_selector)
    ios[1].set_structure(structure_2)
    ios[1].save(interface_path_2, select=structure_2_selector)

    return True

def clean_structure(pdb_path, save_path):

    # Class to select non-HET residues
    class Cleaner(Select):
        def accept_residue(self, residue):
            if one_letter_code(residue.get_resname()) == "":
                return 0
            else:
                return 1

    # Make the right parser and IO objects
    extension = os.path.basename(pdb_path).split(".")[-1]
    if extension == "cif": 
        parser = MMCIFParser(QUIET=True)
        io = MMCIFIO()
    elif extension == "ent": 
        parser = PDBParser(QUIET=True)
        io = PDBIO()

    # Load the structure
    s = parser.get_structure("cleaning", pdb_path)
    # Save the structure
    io.set_structure(s)
    io.save(save_path, select=Cleaner())

    return True
    
        

