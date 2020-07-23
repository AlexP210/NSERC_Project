import My_Library as ml
import os
import sys

if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) < 3:
        print(f"$ Usage: python {os.path.basename(__file__)} <Species Root Directory> <Distance Cutoff>")
        sys.exit()
    root_directory = sys.argv[1]
    distance_cutoff = float(sys.argv[2])

    # Go through each species folder
    for species_folder in os.listdir(root_directory):
        print(species_folder)

        # Set the directories
        species_directory = os.path.join(root_directory, species_folder)
        extracted_chains_directory = os.path.join(species_directory, "Clean_Structures")
        interfaces_directory = os.path.join(species_directory, "Interfaces")
        if not os.path.exists(interfaces_directory): os.mkdir(interfaces_directory)

        # Take each .cif in the Clean_Structures folder, and save the interface to the
        # Interfaces folder

        # Find which chains we have to deal with for each pdb id
        chain_id_map = {}
        print("     Mapping IDs to chain letters ...")
        for filename in os.listdir(extracted_chains_directory):
            structure_name = filename.split(".")[0]
            pdb_id = structure_name.split("_")[0]
            chain_letter = structure_name.split("_")[1]
            if pdb_id not in chain_id_map:
                chain_id_map[pdb_id] = [chain_letter, ]
            else:
                chain_id_map[pdb_id].append(chain_letter)
        
        # Go through each pair of chains, and get the interfaces
        print("     Extracting Interfaces ...")
        for pdb_id, chain_letters in chain_id_map.items():
            for chain_letter_a_idx in range(len(chain_letters)):
                chain_letter_a = chain_letters[chain_letter_a_idx]
                for chain_letter_b in chain_letters[:chain_letter_a_idx]:
                    print(f"         {pdb_id}")
                    # Get the chain and interface file paths
                    chain_a_filename = f"{pdb_id}_{chain_letter_a}.ent"
                    chain_b_filename = f"{pdb_id}_{chain_letter_b}.ent"

                    chain_a_path = os.path.join(extracted_chains_directory, chain_a_filename)
                    chain_b_path = os.path.join(extracted_chains_directory, chain_b_filename)

                    interface_a_filename = f"{pdb_id}{chain_letter_a}{chain_letter_b}_{chain_letter_a}.ent"
                    interface_b_filename = f"{pdb_id}{chain_letter_a}{chain_letter_b}_{chain_letter_b}.ent"
                    interface_a_path = os.path.join(separated_interfaces_directory, interface_a_filename)
                    interface_b_path = os.path.join(separated_interfaces_directory, interface_b_filename)
                    ml.get_interface(chain_a_path, chain_b_path, distance_cutoff, interface_a_path, interface_b_path)


