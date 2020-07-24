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
        complex_interfaces_directory = os.path.join(species_directory, "Complex_Interfaces")
        separated_complex_interfaces_directory = os.path.join(species_directory, "Separated_Complex_Interfaces")
        if not os.path.exists(separated_complex_interfaces_directory): os.mkdir(separated_complex_interfaces_directory)
        for complex_interface_filename in os.listdir(complex_interfaces_directory):
            print( f"       {complex_interface_filename.split('.')[0]}" )
            complex_interface_path = os.path.join(complex_interfaces_directory, complex_interface_filename)
            ml.separate_chains(complex_interface_path, separated_complex_interfaces_directory)

    