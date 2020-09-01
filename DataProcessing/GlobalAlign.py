import My_Library as ml
import os
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: $ python {os.path.basename(__file__)} <PDB Path 1> <PDB Path 2>")
        sys.exit()
    
    pdb_1_path = sys.argv[1]
    pdb_1_name = os.path.basename(pdb_1_path).split(".")[0]
    pdb_2_path = sys.argv[2]
    pdb_2_name = os.path.basename(pdb_2_path).split(".")[0]

    seq_1 = ml.load_sequence(pdb_1_path)
    seq_2 = ml.load_sequence(pdb_2_path)
    
    print()
    nw_score, PID = ml.calculate_percent_identity(pdb_1_name, pdb_1_path, pdb_2_name, pdb_2_path, print_alignment=True)
    print(f"Percent Identity: {PID}")
    print(f"Alignment Score: {nw_score}")
    print(f"Adjusted Alignment Score: {nw_score / max(len(seq_1), len(seq_2))}")
