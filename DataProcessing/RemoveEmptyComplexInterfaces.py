import My_Library as ml
import os
import sys

if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) != 3:
        print(f"Usage: $ python {__file__.split('.')[0]} <Species Root Directory> <Number of Residues>")
        sys.exit()
    root_directory = sys.argv[1]
    n_residues = int(sys.argv[2])

    for species_folder in os.listdir(root_directory):
        print(species_folder)
        species_directory = os.path.join(root_directory, species_folder)
        complex_interfaces_directory = os.path.join(species_directory, "Complex_Interfaces")
        for complex_interface_filename in complex_interfaces_directory:
            complex_interface_path = os.path.join(complex_interfaces_directory, complex_interface_filename)
            try:
                if len(ml.load_sequence(complex_interface_path)) <= n_residues:
                    os.remove(complex_interface_path)
            except:
                os.remove(complex_interface_path)


    

