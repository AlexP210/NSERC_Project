import os
import sys
import pandas as pd

if __name__ == "__main__":
    
    root_directory = sys.argv[1]
    out_file = os.path.join(root_directory, "..", "PercentSignificant.csv")
    with open(out_file, "w") as out:
        out.write("Species,Total_Number,Fraction_Significant\n")
        for species_folder in os.listdir(root_directory):
            print(species_folder)
            species_directory = os.path.join(root_directory, species_folder)
            comparisons_path = os.path.join(species_directory, "ChainsWithSignificances.csv")
            try:
                comparisons = pd.read_csv(comparisons_path)
                pdb_dict = {}
                for index, row in comparisons.iterrows():
                    print(row)
                    if row[10] < 0.05:
                        if row[1] in pdb_dict:
                            pdb_dict[row[1]] += 1
                        else:
                            pdb_dict[row[1]] = 1
                    else:
                        if not row[1] in pdb_dict:
                            pdb_dict[row[1]] = 0
                n_sig = 0
                n_insig = 0
                print(pdb_dict)
                for pdb_id, count in pdb_dict.items():
                    if count == 3:
                        n_sig += 1
                    else:
                        n_insig += 1
                out.write(f"{species_folder},{len(pdb_dict.keys())},{n_sig/(n_sig + n_insig)}\n")
            except pd.errors.EmptyDataError:
                out.write(f"{species_folder},{0},{0}\n")


