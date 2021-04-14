import pandas as pd
import os
import numpy as np


folder_path = r"C:\Users\alexp\2C2"

species_folder = os.path.join(folder_path, "Data")

data = pd.DataFrame()

for species in os.listdir(species_folder):
    print()
    print(species)
    species_directory = os.path.join(species_folder, species)
    significances_path = os.path.join(species_directory, "ChainsWithSignificances2.csv")
    significance_dataframe = pd.read_csv(significances_path)
    print(significance_dataframe.shape)

    species_labels = np.array([species for _ in range(significance_dataframe.shape[0])])
    significance_dataframe["Species"] = species_labels
    print(significance_dataframe)

    data = data.append(significance_dataframe)
    print(data)

data.to_csv(os.path.join(folder_path, "Output.csv"))