# Changes over V1:
# 1. Cleans from Duplicates_Removed

import My_Library as ml
import sys
import os


if __name__ == "__main__":

    # Check arguments:
    if len(sys.argv) != 2:
        print(f"Usage: $ python {os.path.basename(__file__)} <Species Root Directory>")
        sys.exit()
    root_directory = sys.argv[1]

    # Go through each species folder, through each .cif in the Downloaded_Structures folder,
    # and save a cleaned version in the Cleaned_Structures folder
    print("1. CLEANING STRUCTURES")
    for species_folder in os.listdir(root_directory):
        print(species_folder)
        # Set the directories
        species_directory = os.path.join(root_directory, species_folder)
        dirty_structures_directory = os.path.join(species_directory, "Duplicates_Removed")
        clean_structures_directory = os.path.join(species_directory, "Clean_Structures")
        if not os.path.exists(clean_structures_directory): os.mkdir(clean_structures_directory)
        # Go through the dirty ones, and save to clean
        for cif_filename in os.listdir(dirty_structures_directory):
            pdb_id = cif_filename.split(".")[0]
            print(f"     {pdb_id}")
            dirty_cif_path = os.path.join(dirty_structures_directory, cif_filename)
            clean_cif_path = os.path.join(clean_structures_directory, f"{pdb_id}.pdb")
            ml.clean_structure_HOMs(dirty_cif_path, clean_cif_path)
