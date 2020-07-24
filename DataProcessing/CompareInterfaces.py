import My_Library as ml
import os
import sys
import pandas as pd

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) < 3:
        print(f"Usage: $ python {__file__.split('.')[-1]} <Species Root Folder> <Symmetry Groups>")
    root_directory = sys.argv[1]
    symmetry_groups = sys.argv[2:]

    for species_folder in os.listdir(root_directory):
        print(species_folder)
        # Set the directories
        species_directory = os.path.join(root_directory, species_folder)
        interface_comparisons_path = os.path.join(species_directory, "Interfaces.csv")
        complex_interfaces_directory = os.path.join(species_directory, "Complex_Interfaces")
        separated_complex_interfaces_directory = os.path.join(species_directory, "Separated_Complex_Interfaces")
        pairwise_interfaces_directory = os.path.join(species_directory, "Pairwise_Interfaces")
        temp_directory = os.path.join(species_directory, "_Temp")
        # Create the dataframe
        interface_comparisons = pd.DataFrame()
        # Start iterating through the list of complex interfaces
        for complex_interface_filename in os.listdir(complex_interfaces_directory):
            pdb_id, extension = complex_interface_filename.split(".")
            complex_interface_path = os.path.join(complex_interfaces_directory, complex_interface_filename)
            chain_letters = ml.get_chain_ids(complex_interface_path)
            # Get the complex-level symmetry of the interface
            complex_level_symmetries = ml.calculate_symmetry_rmsd(
                complex_interface_path, 
                symmetry_groups, 
                os.path.join(temp_directory, "AnAnaS_Output.txt"))
            # Go through each pair of interfaces in the complex
            for chain_a_idx in range(1, len(chain_letters)):
                for chain_b_idx in range(chain_a_idx):
                    print(f"Processing {pdb_id}_{chain_a_letter} VS {pdb_id}_{chain_b_letter}")
                    # Get the paths for the separated interfaces
                    chain_a_letter = chain_letters[chain_a_idx]
                    chain_b_letter = chain_letters[chain_b_idx]
                    chain_a_interface_filename = f"{pdb_id}_{chain_a_letter}.{extension}"
                    chain_b_interface_filename = f"{pdb_id}_{chain_b_letter}.{extension}"
                    chain_a_interface_path = os.path.join(separated_complex_interfaces_directory, chain_a_interface_filename)
                    chain_b_interface_path = os.path.join(separated_complex_interfaces_directory, chain_a_interface_filename)
                    
                    # Get the features
                    interface_a_length = len( ml.load_sequence(chain_a_interface_path) )
                    interface_b_length = len( ml.load_sequence(chain_b_interface_path) )
                    length_difference = abs( interface_a_length - interface_b_length )
                    alignment_score, percent_identity = ml.calculate_percent_identity(
                        f"{pdb_id}_{chain_a_letter}",
                        chain_a_interface_path, 
                        f"{pdb_id}_{chain_b_letter}", 
                        chain_b_interface_path, 
                        save_alignment_path = os.path.join(temp_directory, "Global_Alignment.fasta"))
                    adjusted_alignment_score = alignment_score / max(interface_a_length, interface_b_length)
                    tm_score = ml.calculate_TMScore(
                        chain_a_interface_path, 
                        chain_b_interface_path, 
                        os.path.join(temp_directory, "TMAlign_Output.txt"), 
                        alignment = os.path.join(temp_directory, "Global_Alignment.fasta"))
                    possible_pairwise_interface_filenames = [ 
                        f"{pdb_id}{chain_a_letter}{chain_b_letter}.{extension}",
                        f"{pdb_id}{chain_b_letter}{chain_a_letter}.{extension}" ]
                    possible_pairwise_interface_paths = [os.path.join(pairwise_interfaces_directory, fn) for fn in possible_pairwise_interface_filenames]
                    for path in possible_pairwise_interface_paths:
                        if os.path.exists(path): pairwise_interface_path = path
                    pairwise_symmetry = ml.calculate_symmetry_rmsd(
                        pairwise_interface_path, 
                        "c2", 
                        os.path.join(temp_directory, "C2_AnAnaS_Output.txt"))
                    row = {
                        "PDB_ID":pdb_id,
                        "Chain_1":chain_a_letter,
                        "Chain_2":chain_b_letter,
                        "Chain_1_Length":interface_a_length,
                        "Chain_2_Length":interface_b_length,
                        "Length_Difference":length_difference,
                        "Alignment_Score":alignment_score,
                        "Alignment_Score_Adjusted":adjusted_alignment_score,
                        "PID":percent_identity,
                        "TM":tm_score,
                        "c2_RMSD":pairwise_symmetry
                    }
                    for symmetry_idx in range(len(symmetry_groups)):
                        row[symmetry_groups[symmetry_idx]] = complex_level_symmetries[symmetry_idx] 
                    print(row)
                    interface_comparisons.append(row)

        interface_comparisons.to_csv(interface_comparisons_path)
                    


