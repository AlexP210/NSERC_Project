import My_Library as ml
import os
import sys

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) != 2:
        print(f"Usage: python {__file__.split('.')[0]} <Species Root Folder>")
        sys.exit()
    root_directory = sys.argv[1]

    for species_folder in os.listdir(root_directory):
        print(species_folder)
        species_directory = os.path.join(root_directory, species_folder)
        pairwise_interfaces_directory = os.path.join(species_directory, "Pairwise_Interfaces")
        separated_pairwise_interfaces_directory = os.path.join(species_directory, "Separated_Pairwise_Interfaces")
        if not os.path.exists(separated_pairwise_interfaces_directory): os.mkdir(separated_pairwise_interfaces_directory)
        for pairwise_interface_filename in os.listdir(pairwise_interfaces_directory):
            print( f"       {pairwise_interface_filename.split('.')[0]}" )
            pairwise_interface_path = os.path.join(pairwise_interfaces_directory, pairwise_interface_filename)
            # May get value error next, because not every pair of chains in every structure has an interface
            try:
                ml.separate_chains(pairwise_interface_path, separated_pairwise_interfaces_directory)
                print( f"       {pairwise_interface_filename.split('.')[0]}" )
            except ValueError:
                pass

    