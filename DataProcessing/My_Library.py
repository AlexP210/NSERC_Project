from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
from Bio.PDB import MMCIFIO, PDBIO, NeighborSearch
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
import pandas as pd
from sklearn.decomposition import PCA


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
    # Temp fix:
    if sequence_a == "" or sequence_b == "":
        return None
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

def load_chain_sequence(structure, chain_letter):
    seq = ""
    for chain in structure.get_chains():
        if chain.get_id() == chain_letter:
            for residue in chain.get_residues():
                seq += one_letter_code(residue.get_resname())
    return seq

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

def calculate_TMScore(pdb_1_path, pdb_2_path, temp_path, alignment = None, matrix_out = None):
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
    if matrix_out != None:
        call += f" -m {matrix_out}"
    call += f" > {temp_path}"
    os.system(call)
    with open(temp_path, "r") as alignment_output:
        TM_scores = []
        for line in alignment_output:
            if line[0:2] == "Al" or line[0:5] == "ERROR":
                # This is the case that the structure alignment failed (chain alignment is too small)
                pass
            if line[0:len("TM-score= ")] == "TM-score= ":
                TM_score_data = line.split("= ")
                TM_score = float(TM_score_data[1].split(" ")[0])
                TM_scores.append(TM_score)
    if len(TM_scores) != 0:
        TM_score = min(TM_scores)
    else:
        TM_score = -1
    # Why do i have this next line:
    # os.remove(temp_path)
    return TM_score

def get_chain_ids(path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(path.split("/")[-1].split("-")[0], path)
    chain_letters = []
    for chain in structure.get_chains():
        chain_letters.append(chain.get_id())
    return chain_letters

def calculate_percent_identity(pdb_1_name, pdb_1_path, pdb_2_name, pdb_2_path, save_alignment_path = None, print_alignment=False):
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
    aligner.substitution_matrix = pam70 # Biopython bug; doesn't initialize the aligner with the subs mat unless we do this
    # Get the sequences
    sequence_a = load_sequence(pdb_1_path)
    sequence_b = load_sequence(pdb_2_path)
    # If the interface is too small, then return NaN's
    if len(sequence_a) == 0 or len(sequence_b) == 0:
        return (np.nan, np.nan)
    # Perform the alignment adn check the identity
    alignment = aligner.align(sequence_a, sequence_b)[0]
    aligned_subsequences = alignment.aligned
    nw_score = alignment.score
    if print_alignment: print(alignment)
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
        if total / len(X) >= 1-prob:
            return i[0]

def create_cdf(X):
    if len(X) > 0:
        return lambda x: sum([i <= x for i in X])/len(X)
    else: 
        return lambda x: 0

def get_interfaces(pdb_path, distance_cutoff, interface_save_path, chain_letters):

    # Selector class
    class SelectInterface(Select):
        def __init__(self, structure, distance_cutoff, chain_letters=None):
            if chain_letters == None: chain_letters = list([list(structure.get_chains())[i].get_id() for i in range(len(list(structure.get_chains())))])
            self.atom_lists = {}
            self.accepted_residues = set()
            self.distance_cutoff = distance_cutoff
            for chain_letter in chain_letters:
                self.atom_lists[chain_letter] = list( [atom for atom in structure.get_atoms() if (atom.get_parent().get_parent().get_id() != chain_letter and atom.get_parent().get_parent().get_id() in chain_letters)])
        def accept_residue(self, residue):
            # If we already cleared this residue, then include it
            # if residue in self.accepted_residues:
            #     return 1
            # otherwise, get the coordinates, and check if anything atoms in the
            # other chains are within range
            chain_letter = residue.get_parent().get_id()
            atom_list = self.atom_lists[chain_letter]
            searcher = NeighborSearch(atom_list)
            # Go through each atom in the residue, find whats 
            for atom in residue.get_atoms():
                position = atom.get_coord()
                contact_residues = searcher.search(position, self.distance_cutoff, level="R")
                for residue in contact_residues:
                    self.accepted_residues.add(residue)
                if len(contact_residues) > 0: 
                    self.accepted_residues.add(residue)
                    return 1
            return 0
    
    # prep the parser and io objects
    extension = os.path.basename(pdb_path).split(".")[-1]
    if extension == "cif":
        parser = MMCIFParser(QUIET=True)
        io = MMCIFIO()
    elif extension in ("pdb", "ent"):
        parser = PDBParser(QUIET=True)
        io = PDBIO()
    
    # Load in the structure, and save the interface
    structure = parser.get_structure("structure", pdb_path)
    selector = SelectInterface(structure, distance_cutoff)
    io.set_structure(structure)
    io.save(interface_save_path, select=selector)

    return True

def get_interface_2files(pdb_path_1, pdb_path_2, distance_cutoff, interface_path_1, interface_path_2):

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

def clean_structure_HOMs(pdb_path, save_path):

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
        io = PDBIO()
    elif extension == "ent": 
        parser = PDBParser(QUIET=True)
        io = PDBIO()

    # Load the structure
    s = parser.get_structure("cleaning", pdb_path)
    # Save the structure
    io.set_structure(s)
    io.save(save_path, select=Cleaner())

    return True
    
def statisticalize(pca, dataframe, features, pre_pca_means, pre_pca_stdevs, post_pca_means, post_pca_stdevs):
    # # Log transform
    # log_col_names = [f"Log({feature})" for feature in features]
    # for i in range(len(log_col_names)):
    #     feature_name = features[i]
    #     log_feature_name = log_col_names[i] 
    #     if feature_name in ("Alignment_Score"): # NOTE: Not logging because it can be negative
    #         dataframe[log_feature_name] = dataframe[feature_name]
    #     elif feature_name in ("Length_Difference"):
    #         dataframe[log_feature_name] = np.log(dataframe[feature_name] + 1)
    #     elif feature_name in ("PID", "TM"): # NOTE: not logging because max changes in the two datasets, gives inconsistent z scores
    #         dataframe[log_feature_name] = dataframe[feature_name]

    # Normalize
    pre_normalized_col_names = [f"{name}_Normalized" for name in features]
    for i in range(len(pre_normalized_col_names)):
        feature_name = features[i]
        normalized_log_col_name = pre_normalized_col_names[i]
        dataframe[normalized_log_col_name] = ( dataframe[feature_name] - pre_pca_means[i] ) / pre_pca_stdevs[i]

    # Do PCA
    transformed_features = pca.transform(dataframe[pre_normalized_col_names])
    transformed_feature_names = [f"PC{i}" for i in range(len(features))]
    for i in range(len(transformed_feature_names)):
        pc_name = transformed_feature_names[i]
        dataframe[pc_name] = transformed_features[:,i]

    # Normalize PC scores
    normalized_pc_names = [f"{pc}_Normalized" for pc in transformed_feature_names]
    for i in range(len(normalized_pc_names)):
        pc = transformed_feature_names[i]
        normalized_pc_name = normalized_pc_names[i]
        dataframe[normalized_pc_name] = ( dataframe[pc] - post_pca_means[i] ) / post_pca_stdevs[i]

    # Z scores
    sum_of_squares = np.zeros((len(dataframe[features[0]])))
    for name in normalized_pc_names:
        sum_of_squares += dataframe[name]**2

    dataframe["Distances"] = np.sqrt( sum_of_squares )
    
    return dataframe

def separate_chains(pdb_file_path, destination_folder):
    extension = pdb_file_path.split(".")[-1]
    pdb_id = os.path.basename(pdb_file_path).split(".")[0]
    if extension == "cif":
        parser = MMCIFParser(QUIET=True)
    elif extension in ("pdb", "ent"):
        parser = PDBParser(QUIET=True)
    # Always use PDBIO() to write separated chains as .ent files. Otherwise TMAlign requires header information
    # (most likely whats causing the issue) either way, using .ent format for the structureal alignment works fine
    io = PDBIO()

    class SplitChains(Select):
        def __init__(self, chain_letter):
            self.chain_letter = chain_letter
        def accept_chain(self, chain):
            if chain.get_id() == self.chain_letter:
                return 1
            else:
                return 0
    
    chain_letters = get_chain_ids(pdb_file_path)
    structure = parser.get_structure(pdb_id, pdb_file_path)
    io.set_structure(structure)
    for chain_letter in chain_letters:
        selector = SplitChains(chain_letter)
        filename = f"{pdb_id}_{chain_letter}.ent"
        save_path = os.path.join(destination_folder, filename)
        io.save(save_path, select=selector)

    return True

def increment_index_state(index_state, species_for_chains):
    index_state = list(index_state)
    for idx_idx in range(len(index_state)):
        if index_state[idx_idx] < len(species_for_chains[idx_idx]) - 1:
            index_state[idx_idx] += 1
            break
        else:
            index_state[idx_idx] = 0
    return tuple(index_state)
 
def find_species(species_for_chains):
    # # # This function is generalized to work for n-mers, where we don't know how many sub-lists there are # # #
    # Can try implementing recursion later

    # Check if we actually have species for all the chains. If not, return False, None to avoid IndexError
    if sum([species_list == [] for species_list in species_for_chains]) > 0: return False, None
    # 1. Pick a species for the first chain, call this species A
    # 2. Look at all the species in the next chain, see if any of them match species A
    # 3. If yes, go to the next chain. If no, go back to the previous chain and pick the next species
    # 4. If we can "find a path" through all the chains, then return True, and that species.
    # 5. If we do not find a path, return False and None
    index_state = tuple( (0 for _ in range(len(species_for_chains))) )
    initial = True
    while initial or index_state != tuple( (0 for _ in range(len(species_for_chains))) ):
        initial = False
        species = [species_for_chains[idx_idx][index_state[idx_idx]] for idx_idx in range(len(index_state))]
        initial_specie = species[0]
        found_species = True
        for specie in species:
            if specie != initial_specie:
                index_state = increment_index_state(index_state, species_for_chains)
                found_species = False
                break
        if found_species:
            return True, initial_specie
    return False, None

if __name__ == "__main__":
    print(find_species([["A", "B"], ["B", "A"]]))
    print(find_species([[],[]]))