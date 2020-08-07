import My_Library as ml
import os
import sys
import pandas as pd

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) < 3:
        print(f"Usage: $ python {__file__.split('.')[-1]} <Species Root Folder> <Symmetry Groups>")
        sys.exit()
    root_directory = sys.argv[1]
    symmetry_groups = sys.argv[2:]

    for species_folder in os.listdir(root_directory):
        print(species_folder)
        # Set the directories
        species_directory = os.path.join(root_directory, species_folder)
        interface_comparisons_path = os.path.join(species_directory, "Interfaces.csv")
        complex_interfaces_directory = os.path.join(species_directory, "Complex_Interfaces")
        separated_complex_interfaces_directory = os.path.join(species_directory, "Separated_Complex_Interfaces")
        temp_directory = os.path.join(species_directory, "_Temp")
        # Create the dataframe
        columns=["PDB_ID","Chain_1","Chain_2","Chain_1_Length","Chain_2_Length","Length_Difference","Alignment_Score","Alignment_Score_Adjusted","PID","TM"]
        for symmetry_group in symmetry_groups:
            columns.append(f"{symmetry_group}_RMSD")
        interface_comparisons = pd.DataFrame(columns=columns)
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
                    # Get the paths for the separated interfaces
                    chain_a_letter = chain_letters[chain_a_idx]
                    chain_b_letter = chain_letters[chain_b_idx]
                    print(f"Processing {pdb_id}_{chain_a_letter} VS {pdb_id}_{chain_b_letter}")
                    chain_a_interface_filename = f"{pdb_id}_{chain_a_letter}.ent"
                    chain_b_interface_filename = f"{pdb_id}_{chain_b_letter}.ent"
                    chain_a_interface_path = os.path.join(separated_complex_interfaces_directory, chain_a_interface_filename)
                    chain_b_interface_path = os.path.join(separated_complex_interfaces_directory, chain_b_interface_filename)
                    
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
                    }
                    for symmetry_idx in range(len(symmetry_groups)):
                        row[f"{symmetry_groups[symmetry_idx]}_RMSD"] = complex_level_symmetries[symmetry_idx] 
                    print(row)
                    interface_comparisons = interface_comparisons.append(row, ignore_index = True)

        interface_comparisons.to_csv(interface_comparisons_path)
                    


