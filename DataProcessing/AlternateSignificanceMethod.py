import os
import sys
import pandas as pd

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: $ python {os.path.basename(__file__)} <Species Root Directory>")
    species_root_directory = sys.argv[1]

    for species_folder in os.listdir(species_root_directory):
        species_directory = os.path.join(species_root_directory, species_folder)
        chains_csv_path = os.path.join(species_directory, "Chains.csv")
        chains_dataframe = pd.from_csv(chains_csv_path)
        
