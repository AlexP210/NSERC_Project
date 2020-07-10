import os
import sys
import shutil

if __name__ == "__main__":

    # Check Arguments
    if len(sys.argv) != 2:
        print(f"Usage: $ python {os.path.basename(__file__)} <Species Root Directory>")
        sys.exit()
    root_directory = sys.argv[1]

    # Go through each species folder, and delete the folders from anything after the download
    for species_folder in os.listdir(root_directory):
        species_directory = os.path.join(root_directory, species_folder)
        for folder in os.listdir(species_directory):
            directory = os.path.join(species_directory, folder)
            if os.path.isdir(directory) and folder != "Structure_Downloads":
                shutil.rmtree(directory)