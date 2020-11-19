
import sys
import os
import shutil

if __name__ == "__main__":
    root_directory = sys.argv[1]
    for species_fol in os.listdir(root_directory):
        species_dir = os.path.join(root_directory, species_fol)
        sanitized = ""
        for c in species_fol:
            if c in"_ABCDEDGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz":
                sanitized += c
        sanitized_dir = os.path.join(root_directory, sanitized)
        shutil.move(species_dir, sanitized_dir)
        # shutil.rmtree(species_dir)
        
