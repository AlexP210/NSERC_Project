from os.path import join
from os import system
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"\nUsage: $ python {__file__.split('/')[-1]} <chain filename 1> <chain filename 2>\n")
        sys.exit()
    path_1 = sys.argv[1]
    path_2 = sys.argv[2]

    chain_a_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/{path_1}"
    chain_b_path = f"../../mnt/e/NSERC_Data/Extracted_Chains/{path_2}"
    output_path = f"../../mnt/e/NSERC_Data/Selected_Alignments/struc_alignment.txt"
    matrix_path = f"../../mnt/e/NSERC_Data/Selected_Alignments/matrix.txt"
    alignment_path = f"../../mnt/e/NSERC_Data/Selected_Alignments/AlignmentPDB"
    tm_align_path = f"./Protein_Complexes/TMAlign/TMalign"
    output_folder = f"../../mnt/e/NSERC_Data/Selected_Alignments"
    call = f'ubuntu run "cd ~;{tm_align_path} {chain_a_path} {chain_b_path} -m {matrix_path} -o {alignment_path} > {output_path}"'
    system(call)