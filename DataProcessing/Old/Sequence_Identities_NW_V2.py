# Changes over V1:
# 1. Calculates percent identity using sequences from original PDB file
#    Allows us to know that we are aligning the right way around

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
import os.path
import os
import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import align, format_alignment

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

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) != 5:
        print(f"Usage: $ python {__file__} <PDB dir> <chains dir> <chains.txt> <output>")
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
        output_file.write("PDB_ID, Chain_1, Chain_2, Percent_Identity_of_Global_Alignment, TM_Score, Identity_Of_Structural_Alignment\n")
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
                        # Calculate structure alignment
                        chain_a_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/pdb{pdb_id}_{chain_letters[sequence_a_idx]}.ent"
                        chain_b_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/pdb{pdb_id}_{chain_letters[sequence_b_idx]}.ent"
                        output_path = f"../../mnt/c/Users/alexp/Desktop/School/3rd_Year/NSERC/StructuralSequenceSimilarity/Scripts/Stuctural_Similarity/struc_alignment.txt"
                        tm_align_path = f"./Protein_Complexes/TMAlign/TMalign"
                        call = f'ubuntu run "cd ~;{tm_align_path} {chain_a_path} {chain_b_path}>{output_path}"'
                        os.system(call)
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
                        os.remove('C:\\Users\\alexp\\Desktop\\School\\3rd_Year\\NSERC\\StructuralSequenceSimilarity\\Scripts\\Stuctural_Similarity\\struc_alignment.txt')
                        TM_score = min(TM_scores)
                        print(f"{pdb_id}, {chain_letters[sequence_a_idx]}, {chain_letters[sequence_b_idx]}, {seq_score}, {TM_score}, {identity}")
                        output_file.write(f"{pdb_id}, {chain_letters[sequence_a_idx]}, {chain_letters[sequence_b_idx]}, {seq_score}, {TM_score}, {identity}\n")

                    


    




