# Changes over V1:
# 1. Does not create separated interfaces, rather the interfaces in the same file

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
        clean_structures_directory = os.path.join(species_directory, "Clean_Structures")
        interfaces_directory = os.path.join(species_directory, "Complex_Interfaces")
        if not os.path.exists(interfaces_directory): os.mkdir(interfaces_directory)

        # Take each .cif in the Clean_Structures folder, and save the interface to the
        # Interfaces folder
        for structure_file_name in os.listdir(clean_structures_directory):
            structure_file_path = os.path.join(clean_structures_directory, structure_file_name)
            interface_file_path = os.path.join(interfaces_directory, structure_file_name)
            if ml.get_interfaces(structure_file_path, distance_cutoff, interface_file_path):
                print(f"     {structure_file_name.split('.')[0]}")


